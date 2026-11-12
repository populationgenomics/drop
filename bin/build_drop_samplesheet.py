#!/usr/bin/env python3

"""
Build the DROP sample sheet for one Metamist RNA cohort, plus the SG map the seqr export needs.

Run by the BUILD_DROP_SAMPLESHEET process in prepare.nf (see docs/CPG_METAMIST_SEQR_DESIGN.md, T2.1):
- the cohort is looked up directly, so it may span projects; missing, inactive and empty cohorts fail the run
- each RNA SG uses its newest completed, active CRAM analysis from cpg-flow-rdrnaseq; SGs without one are
  skipped with a warning, and CRAM locations are used exactly as recorded in Metamist
- RNA_ID is the SG ID; DROP_GROUP is <sample type>__<SG type>, e.g. lcl__polyarna
- a participant appearing twice in one DROP_GROUP fails the run, as seqr can hold only one of them
- groups of min_group_size or fewer samples are warned about: DROP needs more than 30 for outlier calling

Writes into --output-dir:
    samplesheet.tsv   DROP input (assets/schema_input.json)
    sg_map.tsv        SG -> participant, seqr project, seqr tissue / product codes, exome project membership

    build_drop_samplesheet.py --cohort COH0001 --seqr-data-type exome --output-dir .
"""

import csv
import logging
from argparse import ArgumentParser
from collections import defaultdict
from pathlib import Path

from metamist.graphql import gql, query

logger = logging.getLogger(__name__)

SAMPLESHEET_COLUMNS = [
    'RNA_ID',
    'RNA_BAM_FILE',
    'RNA_BAI_FILE',
    'DNA_VCF_FILE',
    'DNA_TBI_FILE',
    'DNA_ID',
    'DROP_GROUP',
    'PAIRED_END',
    'COUNT_MODE',
    'COUNT_OVERLAPS',
    'STRAND',
    'HPO_TERMS',
    'GENE_COUNTS_FILE',
    'GENE_ANNOTATION',
    'GENOME',
    'SPLICE_COUNTS_DIR',
]
SG_MAP_COLUMNS = [
    'sg_id',
    'participant_eid',
    'metamist_project',
    'seqr_project_guid',
    'sample_type',
    'sg_type',
    'seqr_tissue',
    'seqr_product',
    'in_exome_project',
]

RNA_SG_TYPES = {'polyarna', 'totalrna', 'transcriptome'}
# CRAMs registered by cpg-flow-rdrnaseq's TrimAlignRNA; older registrations carry no stage
CRAM_STAGES = {'TrimAlignRNA', None}
# Metamist sample type -> seqr RnaSample.tissue_type; anything else (pbmc, *-chx) isn't exported
SEQR_TISSUES = {'fibroblast': 'F', 'muscle': 'M', 'lcl': 'L'}
# Metamist SG type -> seqr RnaSample.sequencing_type, added on cardinal-seqr dev (T7.1)
SEQR_PRODUCTS = {'polyarna': 'P', 'totalrna': 'R'}

# fixed for every cpg-flow-rdrnaseq CRAM: paired FASTQs only, reverse-stranded libraries, chr-prefixed contigs
STRAND = 'reverse'
PAIRED_END = 'true'
COUNT_MODE = 'IntersectionStrict'
COUNT_OVERLAPS = 'true'
GENOME = 'ucsc'

COHORT_QUERY = """
query DropSamplesheet($cohort: String!) {
    cohorts(id: {eq: $cohort}) {
        id
        status
        sequencingGroups {
            id
            type
            analyses(type: {eq: "cram"}, status: {eq: COMPLETED}, active: {eq: true}) {
                id
                output
                meta
                timestampCompleted
            }
            sample {
                type
                project {
                    name
                    meta
                }
                participant {
                    externalId
                    samples {
                        sequencingGroups {
                            id
                            type
                            technology
                        }
                    }
                }
            }
        }
    }
}
"""


class SamplesheetError(Exception):
    """An unusable cohort, or sequencing groups DROP or seqr can't take."""


def fetch_cohort_sgs(cohort_id: str) -> list[dict]:
    """The cohort's SGs. Fails if the cohort is missing, inactive, or empty."""
    cohorts = query(gql(COHORT_QUERY), variables={'cohort': cohort_id})['cohorts']
    if not cohorts:
        raise SamplesheetError(f'Cohort {cohort_id} not found')
    cohort = cohorts[0]
    if cohort['status'] != 'active':
        raise SamplesheetError(f'Cohort {cohort_id} is {cohort["status"]}, not active')
    if not cohort['sequencingGroups']:
        raise SamplesheetError(f'Cohort {cohort_id} has no sequencing groups')
    return cohort['sequencingGroups']


def latest_cram(sg: dict) -> str | None:
    """The newest completed, active TrimAlignRNA CRAM, or None."""
    crams = [
        a
        for a in sg['analyses']
        if (a.get('output') or '').endswith('.cram') and (a.get('meta') or {}).get('stage') in CRAM_STAGES
    ]
    if not crams:
        return None
    return max(crams, key=lambda a: (a.get('timestampCompleted') or '', a['id']))['output']


def has_exome_sg(participant: dict) -> bool:
    """The seqr exome project is built from short-read exome SGs, so only these participants are in it."""
    return any(
        other['type'] == 'exome' and other['technology'] == 'short-read'
        for sample in participant.get('samples') or []
        for other in sample['sequencingGroups']
    )


def build_rows(sgs: list[dict], seqr_data_type: str, min_group_size: int) -> tuple[list[dict], list[dict]]:
    """Sample sheet and SG map rows, each SG once. Reports every bad SG at once, not one per run."""
    sheet: list[dict] = []
    sg_map: list[dict] = []
    seen: set[str] = set()
    participants_by_group: dict[str, dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))

    for sg in sgs:
        sg_id = sg['id']
        if sg_id in seen:
            continue
        seen.add(sg_id)
        if sg['type'] not in RNA_SG_TYPES:
            logger.warning(f'{sg_id}: skipping {sg["type"]}, not an RNA sequencing group')
            continue
        if not (cram := latest_cram(sg)):
            logger.warning(f'{sg_id}: skipping, no completed and active TrimAlignRNA CRAM')
            continue

        sample = sg['sample']
        participant = sample.get('participant') or {}
        participant_eid = participant.get('externalId') or ''
        group = f'{sample["type"]}__{sg["type"]}'
        if participant_eid:
            participants_by_group[group][participant_eid].append(sg_id)

        sheet.append(
            {
                'RNA_ID': sg_id,
                'RNA_BAM_FILE': cram,
                'RNA_BAI_FILE': f'{cram}.crai',
                'DROP_GROUP': group,
                'PAIRED_END': PAIRED_END,
                'COUNT_MODE': COUNT_MODE,
                'COUNT_OVERLAPS': COUNT_OVERLAPS,
                'STRAND': STRAND,
                'GENOME': GENOME,
            },
        )
        project = sample['project']
        sg_map.append(
            {
                'sg_id': sg_id,
                'participant_eid': participant_eid,
                'metamist_project': project['name'],
                'seqr_project_guid': (project.get('meta') or {}).get(f'seqr-project-{seqr_data_type}') or '',
                'sample_type': sample['type'],
                'sg_type': sg['type'],
                'seqr_tissue': SEQR_TISSUES.get(sample['type'], ''),
                'seqr_product': SEQR_PRODUCTS.get(sg['type'], ''),
                'in_exome_project': str(has_exome_sg(participant)).lower(),
            },
        )

    errors = [
        f'{group}: participant {eid} has several sequencing groups ({", ".join(sorted(sg_ids))})'
        for group, by_participant in sorted(participants_by_group.items())
        for eid, sg_ids in sorted(by_participant.items())
        if len(sg_ids) > 1
    ]
    if errors:
        raise SamplesheetError('\n'.join(errors))
    if not sheet:
        raise SamplesheetError('No RNA sequencing groups with a CRAM')

    warn_about_seqr_mapping(sg_map, seqr_data_type)
    warn_about_group_sizes(sheet, min_group_size)
    return sheet, sg_map


def warn_about_seqr_mapping(sg_map: list[dict], seqr_data_type: str) -> None:
    """These SGs still go through DROP, but the seqr export will skip them."""
    checks = {
        'no participant': lambda row: not row['participant_eid'],
        f'no seqr-project-{seqr_data_type} in project meta': lambda row: not row['seqr_project_guid'],
        'no seqr tissue for its sample type': lambda row: not row['seqr_tissue'],
        'no seqr sequencing product for its SG type': lambda row: not row['seqr_product'],
        'participant has no exome sequencing group': lambda row: row['in_exome_project'] != 'true',
    }
    for reason, check in checks.items():
        if sg_ids := sorted(row['sg_id'] for row in sg_map if check(row)):
            logger.warning(f'Not exportable to seqr ({reason}): {", ".join(sg_ids)}')


def warn_about_group_sizes(sheet: list[dict], min_group_size: int) -> None:
    sizes: dict[str, int] = defaultdict(int)
    for row in sheet:
        sizes[row['DROP_GROUP']] += 1
    for group, size in sorted(sizes.items()):
        if size <= min_group_size:
            logger.warning(f'{group}: {size} samples, outlier calling needs more than {min_group_size}')


def write_tsv(rows: list[dict], columns: list[str], path: Path) -> None:
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, columns, delimiter='\t', restval='', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main(cohort_id: str, seqr_data_type: str, output_dir: str, min_group_size: int) -> None:
    sheet, sg_map = build_rows(fetch_cohort_sgs(cohort_id), seqr_data_type, min_group_size)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(sheet, SAMPLESHEET_COLUMNS, out / 'samplesheet.tsv')
    write_tsv(sg_map, SG_MAP_COLUMNS, out / 'sg_map.tsv')
    groups = sorted({row['DROP_GROUP'] for row in sheet})
    logger.info(f'{cohort_id}: {len(sheet)} SGs in {len(groups)} DROP groups ({", ".join(groups)}) -> {out}')


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--cohort', required=True, help='Metamist cohort ID')
    parser.add_argument(
        '--seqr-data-type',
        default='exome',
        choices=['exome', 'genome'],
        help='which seqr project the RNA attaches to: project meta seqr-project-<type>',
    )
    parser.add_argument('--output-dir', default='.', help='where samplesheet.tsv and sg_map.tsv are written')
    parser.add_argument('--min-group-size', type=int, default=30, help='warn about DROP groups this size or smaller')
    args = parser.parse_args()
    main(args.cohort, args.seqr_data_type, args.output_dir, args.min_group_size)
