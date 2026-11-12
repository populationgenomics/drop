#!/usr/bin/env python3

"""
Turn DROP's results for one DROP group into the files seqr's project RNA-seq upload reads
(docs/CPG_METAMIST_SEQR_DESIGN.md, T4.1; formats in §6.2).

Run by the EXPORT_SEQR process at the end of main.nf, or on its own against published results:
- sample IDs are swapped from Metamist SG IDs to participant external IDs (seqr's Individual_ID), using sg_map.tsv
- rows are split by (Metamist dataset, seqr tissue, seqr sequencing product): each dataset is its own seqr project,
  and seqr takes tissue and product once per upload; SGs without a tissue or product, or whose participant isn't in
  the seqr exome project, are skipped and reported
- each dataset's seqr project GUID comes from its Metamist project meta (seqr-project-<type>, in sg_map.tsv) and is
  listed in export_report.tsv; a dataset without one is still exported, with a warning
- expression outliers keep seqr's five columns; splice outliers are renamed to seqr's column names and FRASER's gene
  symbols are mapped back to gene IDs through DROP's own gene_name_mapping table; TPM keeps its three columns
- gene IDs lose any version suffix; with --seqr-gene-ids, every gene ID must be in seqr's gene table, or the export
  fails and writes nothing, as seqr would reject the whole file

Writes under --output-dir:
    <dataset>/<tissue>-<product>/expression_outliers_-_<cohort>_<arGuid>.tsv.gz
    <dataset>/<tissue>-<product>/splice_outliers_-_<cohort>_<arGuid>.tsv.gz
    <dataset>/<tissue>-<product>/tpm_-_<cohort>_<arGuid>.tsv.gz
    export_report.tsv     files written, with their seqr project GUID and row and sample counts
    export_skipped.tsv    SGs left out, and why

    export_seqr_rna.py --sg-map sg_map.tsv --cohort COH0001 --ar-guid <guid> \\
        --outrider OUTRIDER_results.tsv --fraser results_per_junction.tsv --tpm tpm.tsv.gz \\
        --gene-name-mapping ensembl105.gene_name_mapping.tsv --seqr-gene-ids seqr_gene_ids.txt
"""

import csv
import gzip
import logging
import re
import sys
from argparse import ArgumentParser
from collections import defaultdict
from pathlib import Path

logger = logging.getLogger(__name__)

EXPRESSION = 'expression_outliers'
SPLICE = 'splice_outliers'
TPM = 'tpm'

# seqr's column names (cardinal-seqr seqr/views/utils/dataset_utils.py:48-98)
EXPRESSION_COLUMNS = ['sampleID', 'geneID', 'pValue', 'padjust', 'zScore']
TPM_COLUMNS = ['sample_id', 'gene_id', 'TPM']
# seqr column -> FRASER results_per_junction.tsv column
SPLICE_SOURCE = {
    'chrom': 'seqnames',
    'start': 'start',
    'end': 'end',
    'strand': 'strand',
    'type': 'type',
    'pValue': 'pValue',
    'pAdjust': 'padjust',
    'deltaIntronJaccardIndex': 'deltaPsi',
    'counts': 'counts',
    'meanCounts': 'meanCounts',
    'totalCounts': 'totalCounts',
    'meanTotalCounts': 'meanTotalCounts',
}
SPLICE_COLUMNS = ['sampleID', 'geneID', *SPLICE_SOURCE]

# ENSG00000123456.12, ENSG00000123456.8_5 (GENCODE lift37), ENSG00000123456.5_PAR_Y
GENE_SUFFIX = re.compile(r'(\.\d+)?(_\d+|_PAR_Y)?$')
NO_GENE = {'', 'NA'}


class ExportError(Exception):
    """Inputs seqr would reject, or that can't be exported."""


def open_text(path: Path, mode: str = 'rt'):
    return gzip.open(path, mode, newline='') if str(path).endswith('.gz') else path.open(mode[0], newline='')


def read_rows(path: Path) -> tuple[list[str], list[dict]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        return list(reader.fieldnames or []), list(reader)


def strip_version(gene_id: str) -> str:
    return GENE_SUFFIX.sub('', gene_id) if gene_id.startswith('ENS') else gene_id


def read_sg_map(path: Path) -> dict[str, dict]:
    with path.open(newline='') as handle:
        return {row['sg_id']: row for row in csv.DictReader(handle, delimiter='\t')}


def read_gene_names(path: Path) -> dict[str, str]:
    """DROP's gene_name_mapping table (comma-separated despite its name): unique gene_name -> gene_id."""
    with path.open(newline='') as handle:
        return {row['gene_name']: strip_version(row['gene_id']) for row in csv.DictReader(handle)}


def read_gene_ids(path: Path) -> set[str]:
    with path.open() as handle:
        return {line.strip() for line in handle if line.strip()}


def destination(
    sg_id: str, sg_map: dict[str, dict], ignore_exome_check: bool
) -> tuple[tuple[str, str, str] | None, str]:
    """(dataset, tissue, product) for an SG, or None with the reason it's left out."""
    row = sg_map.get(sg_id)
    if row is None:
        return None, 'not in sg_map.tsv'
    problems = [
        (not row['participant_eid'], 'no participant'),
        (not row['metamist_project'], 'no Metamist dataset'),
        (not row['seqr_tissue'], f'no seqr tissue for sample type {row["sample_type"]}'),
        (not row['seqr_product'], f'no seqr sequencing product for SG type {row["sg_type"]}'),
        (row['in_exome_project'] != 'true' and not ignore_exome_check, 'participant not in the seqr exome project'),
    ]
    for failed, reason in problems:
        if failed:
            return None, reason
    return (row['metamist_project'], row['seqr_tissue'], row['seqr_product']), ''


def seqr_projects(sg_map: dict[str, dict], project_override: str | None) -> dict[str, str]:
    """Metamist dataset -> seqr project GUID (its project meta seqr-project-<type>, carried in sg_map.tsv)."""
    guids: dict[str, set[str]] = defaultdict(set)
    for row in sg_map.values():
        guids[row['metamist_project']].add(project_override or row['seqr_project_guid'])
    if clashes := {dataset: g for dataset, g in guids.items() if len(g - {''}) > 1}:
        raise ExportError(f'datasets with several seqr project GUIDs in sg_map.tsv: {clashes}')
    return {dataset: max(g) for dataset, g in guids.items()}


class Writers:
    """One gzipped TSV per (data type, dataset, tissue, product), opened on first use."""

    def __init__(self, output_dir: Path, cohort: str, ar_guid: str):
        self.output_dir = output_dir
        self.suffix = f'_-_{cohort}_{ar_guid}.tsv.gz'
        self.files: dict[tuple, tuple] = {}
        self.rows: dict[tuple, int] = defaultdict(int)
        self.samples: dict[tuple, set] = defaultdict(set)

    def write(self, data_type: str, columns: list[str], where: tuple[str, str, str], row: dict, sample: str):
        key = (data_type, *where)
        if key not in self.files:
            path = self.output_dir / where[0] / f'{where[1]}-{where[2]}' / f'{data_type}{self.suffix}'
            path.parent.mkdir(parents=True, exist_ok=True)
            # kept open across rows and closed in close(); a context manager can't span the writes
            handle = gzip.open(path, 'wt', newline='')  # noqa: SIM115
            writer = csv.DictWriter(handle, columns, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            self.files[key] = (path, handle, writer)
        self.files[key][2].writerow(row)
        self.rows[key] += 1
        self.samples[key].add(sample)

    def close(self) -> list[Path]:
        for _, handle, _ in self.files.values():
            handle.close()
        return [path for path, _, _ in self.files.values()]


def require(header: list[str], columns: list[str], what: str):
    if missing := [c for c in columns if c not in header]:
        raise ExportError(f'{what}: missing column(s) {", ".join(missing)}')


class Exporter:
    """Collects rows from DROP's tables into seqr files, tracking skipped SGs and gene IDs seqr would reject."""

    def __init__(  # noqa: PLR0913
        self,
        sg_map: dict[str, dict],
        output_dir: Path,
        *,
        cohort: str,
        ar_guid: str,
        gene_names: dict[str, str] | None,
        seqr_gene_ids: set[str] | None,
        project_override: str | None,
        ignore_exome_check: bool,
    ):
        self.sg_map = sg_map
        self.output_dir = output_dir
        self.writers = Writers(output_dir, cohort, ar_guid)
        self.gene_names = gene_names
        self.seqr_gene_ids = seqr_gene_ids
        self.seqr_projects = seqr_projects(sg_map, project_override)
        self.ignore_exome_check = ignore_exome_check
        self.skipped: dict[str, str] = {}
        self.unknown_genes: set[str] = set()
        self.unmapped_symbols: set[str] = set()

    def where_for(self, sg_id: str) -> tuple[str, str, str] | None:
        where, reason = destination(sg_id, self.sg_map, self.ignore_exome_check)
        if where is None:
            self.skipped[sg_id] = reason
        return where

    def check_gene(self, gene_id: str):
        if self.seqr_gene_ids is not None and gene_id not in NO_GENE and gene_id not in self.seqr_gene_ids:
            self.unknown_genes.add(gene_id)

    def participant(self, sg_id: str) -> str:
        return self.sg_map[sg_id]['participant_eid']

    def outrider(self, path: Path):
        header, rows = read_rows(path)
        if rows:
            require(header, EXPRESSION_COLUMNS, path.name)
        for row in rows:
            if (where := self.where_for(row['sampleID'])) is None:
                continue
            gene_id = strip_version(row['geneID'])
            self.check_gene(gene_id)
            participant = self.participant(row['sampleID'])
            out = {c: row[c] for c in EXPRESSION_COLUMNS} | {'sampleID': participant, 'geneID': gene_id}
            self.writers.write(EXPRESSION, EXPRESSION_COLUMNS, where, out, participant)

    def splice_gene_ids(self, symbols: str) -> str:
        """FRASER's ;-joined gene symbols as ;-joined gene IDs (seqr splits on ;), or NA."""
        gene_ids = []
        for raw in (symbols or '').split(';'):
            symbol = raw.strip()
            if symbol in NO_GENE:
                continue
            if symbol not in self.gene_names:
                self.unmapped_symbols.add(symbol)
                continue
            gene_ids.append(self.gene_names[symbol])
            self.check_gene(self.gene_names[symbol])
        return ';'.join(gene_ids) or 'NA'

    def fraser(self, path: Path):
        header, rows = read_rows(path)
        if rows:
            require(header, ['sampleID', 'hgncSymbol', *SPLICE_SOURCE.values()], path.name)
            if self.gene_names is None:
                raise ExportError('--gene-name-mapping is needed to map FRASER gene symbols to gene IDs')
        for row in rows:
            if (where := self.where_for(row['sampleID'])) is None:
                continue
            participant = self.participant(row['sampleID'])
            out = {seqr: row[drop] for seqr, drop in SPLICE_SOURCE.items()}
            out |= {'sampleID': participant, 'geneID': self.splice_gene_ids(row['hgncSymbol'])}
            self.writers.write(SPLICE, SPLICE_COLUMNS, where, out, participant)

    def tpm(self, path: Path):
        header, rows = read_rows(path)
        if rows:
            require(header, TPM_COLUMNS, path.name)
        for row in rows:
            if (where := self.where_for(row['sample_id'])) is None:
                continue
            gene_id = strip_version(row['gene_id'])
            self.check_gene(gene_id)
            participant = self.participant(row['sample_id'])
            out = {'sample_id': participant, 'gene_id': gene_id, 'TPM': row['TPM']}
            self.writers.write(TPM, TPM_COLUMNS, where, out, participant)

    def finish(self) -> list[dict]:
        """Closes the files; on any gene seqr would reject, deletes them and raises."""
        paths = self.writers.close()
        errors = []
        if self.unmapped_symbols:
            errors.append(
                f'FRASER gene symbols not in gene_name_mapping: {", ".join(sorted(self.unmapped_symbols)[:50])}'
            )
        if self.unknown_genes:
            errors.append(
                f"{len(self.unknown_genes)} gene IDs not in seqr's gene table, so seqr would reject the file: "
                f'{", ".join(sorted(self.unknown_genes)[:50])}'
            )
        if errors:
            for path in paths:
                path.unlink(missing_ok=True)
            raise ExportError('\n'.join(errors))
        return [
            {
                'data_type': key[0],
                'dataset': key[1],
                'seqr_project_guid': self.seqr_projects[key[1]],
                'seqr_tissue': key[2],
                'seqr_product': key[3],
                'file': str(path.relative_to(self.output_dir)),
                'rows': self.writers.rows[key],
                'samples': len(self.writers.samples[key]),
            }
            for key, (path, _, _) in self.writers.files.items()
        ]


def export(  # noqa: PLR0913
    sg_map: dict[str, dict],
    cohort: str,
    ar_guid: str,
    output_dir: Path,
    *,
    outrider: Path | None = None,
    fraser: Path | None = None,
    tpm: Path | None = None,
    gene_names: dict[str, str] | None = None,
    seqr_gene_ids: set[str] | None = None,
    project_override: str | None = None,
    ignore_exome_check: bool = False,
) -> tuple[list[dict], dict[str, str]]:
    """Writes the seqr files; returns the report rows and the skipped SGs with reasons."""
    exporter = Exporter(
        sg_map,
        output_dir,
        cohort=cohort,
        ar_guid=ar_guid,
        gene_names=gene_names,
        seqr_gene_ids=seqr_gene_ids,
        project_override=project_override,
        ignore_exome_check=ignore_exome_check,
    )
    if outrider:
        exporter.outrider(outrider)
    if fraser:
        exporter.fraser(fraser)
    if tpm:
        exporter.tpm(tpm)
    return exporter.finish(), exporter.skipped


def write_report(report: list[dict], skipped: dict[str, str], output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    columns = ['data_type', 'dataset', 'seqr_project_guid', 'seqr_tissue', 'seqr_product', 'file', 'rows', 'samples']
    with (output_dir / 'export_report.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, columns, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(report)
    with (output_dir / 'export_skipped.tsv').open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['sg_id', 'reason'])
        writer.writerows(sorted(skipped.items()))


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--sg-map', required=True, type=Path, help='sg_map.tsv from prepare.nf')
    parser.add_argument('--cohort', required=True, help='Metamist cohort ID, used in file names')
    parser.add_argument('--ar-guid', default='local', help='analysis-runner GUID, used in file names')
    parser.add_argument('--outrider', type=Path, help='OUTRIDER_results.tsv')
    parser.add_argument('--fraser', type=Path, help='FRASER results_per_junction.tsv')
    parser.add_argument('--tpm', type=Path, help='tpm.tsv.gz')
    parser.add_argument('--gene-name-mapping', type=Path, help="DROP's <annotation>.gene_name_mapping.tsv")
    parser.add_argument('--seqr-gene-ids', type=Path, help="seqr's gene IDs, one per line; skips the check if absent")
    parser.add_argument('--seqr-project-override', help='report this seqr project GUID for every dataset (testing)')
    parser.add_argument(
        '--ignore-exome-check', action='store_true', help='keep SGs whose participant has no exome SG (testing)'
    )
    parser.add_argument('--output-dir', type=Path, default=Path('seqr_export'))
    args = parser.parse_args()

    if args.seqr_gene_ids is None:
        logger.warning('No --seqr-gene-ids: gene IDs are not checked against seqr')
    try:
        report, skipped = export(
            read_sg_map(args.sg_map),
            args.cohort,
            args.ar_guid,
            args.output_dir,
            outrider=args.outrider,
            fraser=args.fraser,
            tpm=args.tpm,
            gene_names=read_gene_names(args.gene_name_mapping) if args.gene_name_mapping else None,
            seqr_gene_ids=read_gene_ids(args.seqr_gene_ids) if args.seqr_gene_ids else None,
            project_override=args.seqr_project_override,
            ignore_exome_check=args.ignore_exome_check,
        )
    except ExportError as error:
        logger.error(str(error))
        sys.exit(1)
    write_report(report, skipped, args.output_dir)
    for line in report:
        logger.info(
            f'{line["file"]}: {line["rows"]} rows, {line["samples"]} samples, '
            f'seqr project {line["seqr_project_guid"] or "unknown"}'
        )
    for dataset in sorted({line['dataset'] for line in report if not line['seqr_project_guid']}):
        logger.warning(
            f'{dataset}: no seqr-project-<type> in its Metamist project meta, so no seqr project to upload to'
        )
    for sg_id, reason in sorted(skipped.items()):
        logger.warning(f'Skipped {sg_id}: {reason}')
    if not report:
        logger.warning('Nothing written: no rows, or every SG was skipped')


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')
    main()
