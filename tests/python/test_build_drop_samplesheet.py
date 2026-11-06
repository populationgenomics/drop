"""bin/build_drop_samplesheet.py: Metamist responses -> DROP sample sheet and SG map. The query is validated
against the schema bundled with the metamist client; responses are hand-built dicts in the shape it returns."""

import csv
import json
import logging
import re
from pathlib import Path

import build_drop_samplesheet
import pytest
from build_drop_samplesheet import (
    SAMPLESHEET_COLUMNS,
    SG_MAP_COLUMNS,
    SamplesheetError,
    build_rows,
    fetch_cohort_sgs,
    main,
)
from conftest import cohort, cram, other_sg, sg
from metamist.graphql import gql, validate

SCHEMA_INPUT = Path(__file__).parents[2] / 'assets' / 'schema_input.json'


def rows(sgs: list[dict], min_group_size: int = 0) -> tuple[list[dict], list[dict]]:
    return build_rows(sgs, 'exome', min_group_size)


def by_sg(sg_map: list[dict]) -> dict[str, dict]:
    return {row['sg_id']: row for row in sg_map}


def mock_query(monkeypatch: pytest.MonkeyPatch, cohorts: list[dict]) -> None:
    monkeypatch.setattr(build_drop_samplesheet, 'query', lambda *_args, **_kwargs: {'cohorts': cohorts})


def test_query_matches_metamist_schema():
    assert validate(gql(build_drop_samplesheet.COHORT_QUERY), use_local_schema=True)


def test_sheet_rows():
    sheet, _ = rows([sg('CPG1')])
    assert sheet == [
        {
            'RNA_ID': 'CPG1',
            'RNA_BAM_FILE': 'gs://cpg-rdnow-main/transcriptome/cram/CPG1.cram',
            'RNA_BAI_FILE': 'gs://cpg-rdnow-main/transcriptome/cram/CPG1.cram.crai',
            'DROP_GROUP': 'lcl__polyarna',
            'PAIRED_END': 'true',
            'COUNT_MODE': 'IntersectionStrict',
            'COUNT_OVERLAPS': 'true',
            'STRAND': 'reverse',
            'GENOME': 'ucsc',
        },
    ]


def test_sheet_values_pass_drop_input_schema():
    """The patterns and enums DROP's nf-schema validation applies to every row."""
    properties = json.loads(SCHEMA_INPUT.read_text())['items']['properties']
    sheet, _ = rows([sg('CPG1'), sg('CPG2', sample_type='fibroblast', sg_type='totalrna')])
    for row in sheet:
        for column, value in row.items():
            spec = properties[column]
            if 'pattern' in spec:
                assert re.search(spec['pattern'], value), (column, value)
            if 'enum' in spec:
                assert value in spec['enum'], (column, value)


def test_sg_map_rows():
    _, sg_map = rows([sg('CPG1', participant='P1')])
    assert sg_map == [
        {
            'sg_id': 'CPG1',
            'participant_eid': 'P1',
            'metamist_project': 'rdnow',
            'seqr_project_guid': 'R_rdnow_exome',
            'sample_type': 'lcl',
            'sg_type': 'polyarna',
            'seqr_tissue': 'L',
            'seqr_product': 'P',
            'in_exome_project': 'true',
        },
    ]


def test_drop_group_is_sample_type_and_sg_type():
    sheet, _ = rows(
        [
            sg('CPG1', sample_type='fibroblast', sg_type='polyarna'),
            sg('CPG2', sample_type='fibroblast-chx', sg_type='polyarna'),
            sg('CPG3', sample_type='muscle', sg_type='totalrna'),
        ],
    )
    assert [r['DROP_GROUP'] for r in sheet] == [
        'fibroblast__polyarna',
        'fibroblast-chx__polyarna',
        'muscle__totalrna',
    ]


def test_cram_location_is_used_verbatim():
    location = 'gs://cpg-rdnow-test/transcriptome/cram/CPG737320.cram'
    sheet, _ = rows([sg('CPG737320', analyses=[cram(location)])])
    assert sheet[0]['RNA_BAM_FILE'] == location


def test_newest_cram_wins():
    analyses = [
        cram('gs://b/old.cram', analysis_id=9, completed='2025-01-01'),
        cram('gs://b/new.cram', analysis_id=3, completed='2026-06-01'),
    ]
    assert rows([sg('CPG1', analyses=analyses)])[0][0]['RNA_BAM_FILE'] == 'gs://b/new.cram'


def test_crams_from_other_stages_are_ignored():
    analyses = [cram('gs://b/other.cram', stage='SomethingElse'), cram('gs://b/legacy.cram', stage=None)]
    assert rows([sg('CPG1', analyses=analyses)])[0][0]['RNA_BAM_FILE'] == 'gs://b/legacy.cram'


def test_sg_without_cram_is_skipped(caplog):
    with caplog.at_level(logging.WARNING):
        sheet, sg_map = rows([sg('CPG1'), sg('CPG2', analyses=[]), sg('CPG3', analyses=[cram('gs://b/x.bam')])])
    assert [r['RNA_ID'] for r in sheet] == ['CPG1']
    assert [r['sg_id'] for r in sg_map] == ['CPG1']
    assert 'CPG2: skipping, no completed and active TrimAlignRNA CRAM' in caplog.text
    assert 'CPG3: skipping' in caplog.text


def test_non_rna_sg_is_skipped(caplog):
    with caplog.at_level(logging.WARNING):
        sheet, _ = rows([sg('CPG1'), sg('CPG2', sg_type='genome')])
    assert [r['RNA_ID'] for r in sheet] == ['CPG1']
    assert 'CPG2: skipping genome' in caplog.text


def test_sgs_are_deduplicated():
    sheet, sg_map = rows([sg('CPG1'), sg('CPG1')])
    assert len(sheet) == len(sg_map) == 1


def test_participant_twice_in_one_group_fails():
    sgs = [sg('CPG1', participant='P1'), sg('CPG2', participant='P1'), sg('CPG3', participant='P2')]
    with pytest.raises(
        SamplesheetError, match=r'lcl__polyarna: participant P1 has several sequencing groups \(CPG1, CPG2\)'
    ):
        rows(sgs)


def test_participant_in_two_groups_is_fine():
    sheet, _ = rows([sg('CPG1', participant='P1'), sg('CPG2', participant='P1', sample_type='fibroblast')])
    assert len(sheet) == 2


def test_no_usable_sg_fails():
    with pytest.raises(SamplesheetError, match='No RNA sequencing groups with a CRAM'):
        rows([sg('CPG1', analyses=[])])


def test_unmapped_tissue_and_product_stay_in_drop_but_are_flagged(caplog):
    with caplog.at_level(logging.WARNING):
        sheet, sg_map = rows([sg('CPG1', sample_type='pbmc'), sg('CPG2', sg_type='transcriptome')])
    assert len(sheet) == 2
    assert by_sg(sg_map)['CPG1']['seqr_tissue'] == ''
    assert by_sg(sg_map)['CPG2']['seqr_product'] == ''
    assert 'no seqr tissue for its sample type): CPG1' in caplog.text
    assert 'no seqr sequencing product for its SG type): CPG2' in caplog.text


def test_exome_project_membership():
    _, sg_map = rows(
        [
            sg('CPG1'),
            sg('CPG2', participant_sgs=[other_sg('CPG20', sg_type='genome')]),
            sg('CPG3', participant_sgs=[other_sg('CPG30', technology='long-read')]),
            sg('CPG4', participant_sgs=[]),
        ],
    )
    assert {k: v['in_exome_project'] for k, v in by_sg(sg_map).items()} == {
        'CPG1': 'true',
        'CPG2': 'false',
        'CPG3': 'false',
        'CPG4': 'false',
    }


def test_seqr_project_comes_from_each_sgs_own_project():
    _, sg_map = rows([sg('CPG1', project='rdnow'), sg('CPG2', project='mito-mdt'), sg('CPG3', project_meta={})])
    assert {k: v['seqr_project_guid'] for k, v in by_sg(sg_map).items()} == {
        'CPG1': 'R_rdnow_exome',
        'CPG2': 'R_mito-mdt_exome',
        'CPG3': '',
    }


def test_seqr_data_type_picks_the_project_meta_key():
    meta = {'seqr-project-exome': 'R_exome', 'seqr-project-genome': 'R_genome'}
    _, sg_map = build_rows([sg('CPG1', project_meta=meta)], 'genome', 0)
    assert sg_map[0]['seqr_project_guid'] == 'R_genome'


def test_sample_without_participant_is_kept_but_flagged(caplog):
    sgs = [sg('CPG1')]
    sgs[0]['sample']['participant'] = None
    with caplog.at_level(logging.WARNING):
        sheet, sg_map = rows(sgs)
    assert len(sheet) == 1
    assert sg_map[0]['participant_eid'] == ''
    assert sg_map[0]['in_exome_project'] == 'false'
    assert 'Not exportable to seqr (no participant): CPG1' in caplog.text


def test_small_groups_are_warned_about(caplog):
    with caplog.at_level(logging.WARNING):
        rows([sg('CPG1'), sg('CPG2'), sg('CPG3', sample_type='muscle')], min_group_size=1)
    assert 'lcl__polyarna' not in caplog.text
    assert 'muscle__polyarna: 1 samples, outlier calling needs more than 1' in caplog.text


@pytest.mark.parametrize(
    ('cohorts', 'message'),
    [
        ([], 'Cohort COH1 not found'),
        ([cohort([sg('CPG1')], status='archived')], 'Cohort COH1 is archived, not active'),
        ([cohort([])], 'Cohort COH1 has no sequencing groups'),
    ],
)
def test_unusable_cohort_fails(monkeypatch, cohorts, message):
    mock_query(monkeypatch, cohorts)
    with pytest.raises(SamplesheetError, match=message):
        fetch_cohort_sgs('COH1')


def test_main_writes_both_files(monkeypatch, tmp_path):
    mock_query(monkeypatch, [cohort([sg('CPG1', participant='P1'), sg('CPG2', participant='P2')])])
    main('COH1', 'exome', str(tmp_path / 'out'), 30)

    with (tmp_path / 'out' / 'samplesheet.tsv').open() as handle:
        sheet = list(csv.DictReader(handle, delimiter='\t'))
    with (tmp_path / 'out' / 'sg_map.tsv').open() as handle:
        sg_map = list(csv.DictReader(handle, delimiter='\t'))

    assert list(sheet[0].keys()) == SAMPLESHEET_COLUMNS
    assert [r['RNA_ID'] for r in sheet] == ['CPG1', 'CPG2']
    assert sheet[0]['DNA_ID'] == ''
    assert list(sg_map[0].keys()) == SG_MAP_COLUMNS
    assert [r['participant_eid'] for r in sg_map] == ['P1', 'P2']
