"""Metamist-shaped fixtures for the sample sheet tests. Metamist is never contacted."""


def cram(output: str, analysis_id: int = 1, stage: str | None = 'TrimAlignRNA', completed: str = '2026-01-01') -> dict:
    return {'id': analysis_id, 'output': output, 'meta': {'stage': stage}, 'timestampCompleted': completed}


def other_sg(sg_id: str, sg_type: str = 'exome', technology: str = 'short-read') -> dict:
    return {'id': sg_id, 'type': sg_type, 'technology': technology}


def sg(  # noqa: PLR0913
    sg_id: str,
    *,
    participant: str | None = None,
    sample_type: str = 'lcl',
    sg_type: str = 'polyarna',
    project: str = 'rdnow',
    project_meta: dict | None = None,
    analyses: list[dict] | None = None,
    participant_sgs: list[dict] | None = None,
) -> dict:
    """An RNA sequencing group, shaped like the response to COHORT_QUERY. Defaults to an exportable lcl SG."""
    if analyses is None:
        analyses = [cram(f'gs://cpg-{project}-main/transcriptome/cram/{sg_id}.cram')]
    if participant_sgs is None:
        participant_sgs = [other_sg(f'{sg_id}-exome')]
    return {
        'id': sg_id,
        'type': sg_type,
        'analyses': analyses,
        'sample': {
            'type': sample_type,
            'project': {
                'name': project,
                'meta': {'seqr-project-exome': f'R_{project}_exome'} if project_meta is None else project_meta,
            },
            'participant': {
                'externalId': participant or f'P-{sg_id}',
                'samples': [{'sequencingGroups': participant_sgs}],
            },
        },
    }


def cohort(sgs: list[dict], cohort_id: str = 'COH1', status: str = 'active') -> dict:
    return {'id': cohort_id, 'status': status, 'sequencingGroups': sgs}
