#!/bin/bash -l
# ============================================================================
# prepare_rerun.sh - park the previous run's outputs before rebuilding.
#
# The pipeline scripts OVERWRITE their outputs but never DELETE anything, so a
# rerun over a populated tree leaves behind any file the new run does not
# happen to rewrite. That is harmless when the sample set and the stage list
# are identical and every stage completes - and quietly wrong otherwise.
#
# This run is the "otherwise" case: the reference changed, so every INTERMEDIATE
# downstream of Stage 0a is stale by definition.
#
# stage3/ is deliberately NOT touched. Those are results, not intermediates:
# each filter setting writes under its own label and the labels are meant to
# accumulate across runs. Moving them would fight that workflow. Provenance is
# handled instead - Stage 3 now stamps reference_id and schema_version into
# summary.tsv, overview.tsv, overview_long.tsv and params.json, so an old file
# and a new one are always distinguishable even though they sit side by side.
# (Pass --include-stage3 if you do want them parked anyway.)
#
# Archiving the intermediates is OPTIONAL too - Stage 3 verifies that every
# sample's Stage 2a output carries the same reference_id and aborts on a
# mixture, so a half-rerun cannot silently blend two references. The reasons to
# archive anyway are that a `mv` within the same filesystem is instant, and that
# it makes a failed snpEff rebuild loud (see the note at the bottom).
#
# Usage:
#   bash prepare_rerun.sh                     # show what would move (dry run)
#   bash prepare_rerun.sh --go                # move intermediates only
#   bash prepare_rerun.sh --go --include-stage3
#   bash prepare_rerun.sh --go --tag pre_refswitch
# ============================================================================
set -euo pipefail

BASE_DIR="/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private"
TE="${BASE_DIR}/transcription_errors"
CANON="${BASE_DIR}/refgenome/canonical"
SNPEFF_DATA="${BASE_DIR}/refgenome/snpeff_data"
SNPEFF_GENOME="BDGP6.46.canonical"

GO=0
INCLUDE_STAGE3=0
TAG="pre_refswitch_$(date +%Y%m%d)"
while [ $# -gt 0 ]; do
    case "$1" in
        --go)              GO=1 ;;
        --include-stage3)  INCLUDE_STAGE3=1 ;;
        --tag)             shift; TAG="$1" ;;
        *) echo "unknown argument: $1"; exit 1 ;;
    esac
    shift
done

ARCHIVE="${TE}/_archive/${TAG}"

# Intermediates only. Ordered as the pipeline produces them.
TARGETS=(
    "${CANON}"
    "${TE}/stage1"
    "${TE}/stage2a"
    "${TE}/stage2b"
)
if [ "${INCLUDE_STAGE3}" = 1 ]; then
    TARGETS+=("${TE}/stage3")
fi

echo "Archive destination : ${ARCHIVE}"
echo "Mode                : $([ "${GO}" = 1 ] && echo 'MOVE (--go given)' || echo 'dry run')"
echo "stage3/ results     : $([ "${INCLUDE_STAGE3}" = 1 ] && echo 'INCLUDED (--include-stage3)' || echo 'left in place (labelled results accumulate there)')"
echo ""

for d in "${TARGETS[@]}"; do
    if [ -e "${d}" ]; then
        n=$(find "${d}" -type f 2>/dev/null | wc -l)
        sz=$(du -sh "${d}" 2>/dev/null | cut -f1)
        printf '  %-60s %8s files  %6s\n' "${d}" "${n}" "${sz}"
    else
        printf '  %-60s %s\n' "${d}" "(absent)"
    fi
done

# The snpEff database is a special case. Stage 0b rebuilds it under the SAME
# genome name, so a successful rebuild overwrites snpEffectPredictor.bin. But
# if the build FAILS, the previous .bin survives and Stage 2b will happily
# annotate against the old reference without complaining. Moving it forces the
# failure to be loud.
BIN="${SNPEFF_DATA}/${SNPEFF_GENOME}"
if [ -e "${BIN}" ]; then
    printf '  %-60s %s\n' "${BIN}" "(snpEff DB - see note below)"
fi

echo ""
if [ "${GO}" != 1 ]; then
    echo "Dry run. Nothing moved. Re-run with --go to move these into ${ARCHIVE}."
    exit 0
fi

mkdir -p "${ARCHIVE}"
for d in "${TARGETS[@]}"; do
    [ -e "${d}" ] || continue
    dest="${ARCHIVE}/$(basename "${d}")"
    if [ -e "${dest}" ]; then
        echo "SKIP (destination exists): ${dest}"
        continue
    fi
    echo "moving ${d} -> ${dest}"
    mv "${d}" "${dest}"
done

if [ -e "${BIN}" ]; then
    dest="${ARCHIVE}/snpeff_${SNPEFF_GENOME}"
    if [ -e "${dest}" ]; then
        echo "SKIP (destination exists): ${dest}"
    else
        echo "moving ${BIN} -> ${dest}"
        mv "${BIN}" "${dest}"
    fi
fi

echo ""
echo "Done. The old intermediates are at ${ARCHIVE}."
if [ "${INCLUDE_STAGE3}" != 1 ]; then
    echo ""
    echo "stage3/ was left untouched, so your labelled results keep accumulating."
    echo "New files there will carry reference_id and schema_version columns;"
    echo "files written before this upgrade have neither, which is how you tell"
    echo "them apart. Remember the unsuffixed Stage 3 columns are nuclear-only"
    echo "from schema_version 2 onward, so do not compare them across versions."
fi
echo ""
echo "Now run:  ./submit_pipeline.sh samples.txt"
