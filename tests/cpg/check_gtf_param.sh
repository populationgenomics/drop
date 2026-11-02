#!/usr/bin/env bash
# CPG: check --gtf (T3.1) builds exactly the annotation input --gene_annotation does, and that
# gene_name_mapping.tsv is published (T3.2). See docs/CPG_METAMIST_SEQR_DESIGN.md.
#
# Re-runs the test profile with -resume, passing the demo GTF directly under the same name the YAML gives
# it. If --gtf is right, every annotation-dependent task comes from cache. Reuse the T0.4 launch directory
# so the cache is there; without it this is a full test run (about 40 minutes).
#
#   tests/cpg/check_gtf_param.sh [launch_dir]     # default /tmp/drop_t04

set -euo pipefail

PIPELINE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
LAUNCH_DIR="${1:-/tmp/drop_t04}"
mkdir -p "${LAUNCH_DIR}"
cd "${LAUNCH_DIR}"

# DROP's images are amd64-only; Apple Silicon runs them under emulation
PROFILE='test,docker'
if [[ "$(uname -m)" == 'arm64' ]]; then
    PROFILE+=',emulate_amd64'
fi

# fit the test profile into Docker Desktop's default 8 GB VM
cat > local.config <<'EOF'
process.resourceLimits = [cpus: 4, memory: 7.GB, time: 2.h]
EOF

# the same GTF and name as assets/gene_annotation.yaml; the test profile's YAML is switched off
cat > gtf_params.yaml <<'EOF'
gene_annotation: null
gtf: https://github.com/gagneurlab/drop_demo_data/raw/refs/heads/main/Data/gencode_annotation_trunc.gtf
gtf_name: v29
EOF

NXF_SYNTAX_PARSER=v2 nextflow run "${PIPELINE}" \
    -profile "${PROFILE}" -c local.config -params-file gtf_params.yaml \
    --outdir results -ansi-log false -resume | tee run_gtf.log

# with -ansi-log false and -resume, only executed tasks get a [PROCESS ...] line; cached ones aren't listed
echo
echo '== tasks executed (not from cache):'
grep -E '^\[PROCESS ' run_gtf.log | awk '{print $3}' | sed 's/.*://' | sort | uniq -c || echo '   none'

status=0
if grep -E '^\[PROCESS ' run_gtf.log | grep -qE 'PREPROCESSGENEANNOTATION|GENECOUNTS_COUNTREADS|FRASER_ANNOTATEGENES'; then
    echo 'FAIL: annotation-dependent tasks re-ran, so --gtf built a different input than --gene_annotation'
    status=1
else
    echo 'OK: annotation-dependent tasks all came from cache'
fi

if ls results/processed_data/annotation/*.gene_name_mapping.tsv > /dev/null 2>&1; then
    echo "OK: published $(ls results/processed_data/annotation/*.gene_name_mapping.tsv)"
else
    echo 'FAIL: gene_name_mapping.tsv not published to results/processed_data/annotation/'
    status=1
fi
exit "${status}"