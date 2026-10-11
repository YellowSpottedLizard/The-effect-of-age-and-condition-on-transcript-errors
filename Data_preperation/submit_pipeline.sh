#!/bin/bash -l
# ============================================================================
# submit_pipeline.sh - submit the whole transcription-error pipeline at once.
# Uses SLURM --dependency=afterok so each step starts only after its inputs
# finish successfully. 
#
# Dependency graph:
#   0A ref --+--> 0B snpEff DB ---------------\
#            +--> 0C N/S table ---------------+--> 2B annotate --\
#            +--> 1 SNVs (array) -------------/                   +--> 3 (array)
#            +--> 2A depth (array, also after 0C) --> 2A reduce --/
#
# Usage:
#   sh submit_pipeline.sh samples.txt                # run everything
#   sh submit_pipeline.sh samples.txt --from 2B      # resume from a stage
# Stages: 0A 0B 0C 1 2A 2AR 2B 3
# ============================================================================
set -euo pipefail
SAMPLES="${1:?usage: ./submit_pipeline.sh samples.txt [--from STAGE]}"
FROM="0A"
if [ $# -ge 2 ]; then
    [ "$2" = "--from" ] && [ -n "${3:-}" ] || { echo "usage: ./submit_pipeline.sh samples.txt [--from STAGE]"; exit 1; }
    FROM="$3"
fi
N=$(wc -l < "${SAMPLES}")

DIR="/cfs/klemming/projects/snic/snic2021-23-688/nobackup/private/Friberg_2025_TranscriptError_eriks_copy/TE"

# number of SETTINGS entries in JOB_stage3.sh (one array task each)
N_FILTER_SETTINGS=$(awk '/^SETTINGS=\(/{f=1; next} f && /^\)/{exit} f && /^[[:space:]]*"/{n++} END{print n+0}' "${DIR}/JOB_stage3.sh")
[ "${N_FILTER_SETTINGS}" -gt 0 ] || { echo "ERROR: no SETTINGS entries found in JOB_stage3.sh"; exit 1; }

echo "Samples: ${N}   Filter settings: ${N_FILTER_SETTINGS}   Starting from: ${FROM}"

order="0A 0B 0C 1 2A 2AR 2B 3"
case " ${order} " in *" ${FROM} "*) ;; *) echo "ERROR: unknown stage '${FROM}' (valid: ${order})"; exit 1;; esac
active=0; declare -A RUN
for s in $order; do [ "$s" = "$FROM" ] && active=1; RUN[$s]=$active; done

J0A=""; J0B=""; J0C=""; J1=""; J2A=""; J2AR=""; J2B=""
dep() { [ -n "$1" ] && echo "--dependency=afterok:$1" || echo ""; }
# afterok on every non-empty job id given
deps() { local d=""; for j in "$@"; do [ -n "$j" ] && d="${d}:${j}"; done
         [ -n "$d" ] && echo "--dependency=afterok${d}" || echo ""; }

if [ "${RUN[0A]}" = 1 ]; then
    J0A=$(sbatch --parsable "${DIR}/JOB_stage0A.sh");                              echo "0A  ref prep        : ${J0A}"
fi
if [ "${RUN[0B]}" = 1 ]; then
    J0B=$(sbatch --parsable $(dep "${J0A}") "${DIR}/JOB_stage0B.sh");              echo "0B  snpEff DB       : ${J0B}"
fi
if [ "${RUN[0C]}" = 1 ]; then
    J0C=$(sbatch --parsable $(dep "${J0A}") "${DIR}/JOB_stage0C.sh");              echo "0C  N/S table       : ${J0C}"
fi
if [ "${RUN[1]}" = 1 ]; then
    J1=$(sbatch --parsable $(dep "${J0A}") --array=1-"${N}"%64 "${DIR}/JOB_stage1.sh" "${SAMPLES}")
    echo "1   extract (array) : ${J1}"
fi
if [ "${RUN[2A]}" = 1 ]; then
    J2A=$(sbatch --parsable $(deps "${J0A}" "${J0C}") --array=1-"${N}"%64 "${DIR}/JOB_stage2A.sh" "${SAMPLES}")
    echo "2A  depth (array)   : ${J2A}"
fi
if [ "${RUN[2AR]}" = 1 ]; then
    J2AR=$(sbatch --parsable $(dep "${J2A}") "${DIR}/JOB_stage2A_reduce.sh" "${SAMPLES}")
    echo "2A  reduce          : ${J2AR}"
fi
if [ "${RUN[2B]}" = 1 ]; then
    J2B=$(sbatch --parsable $(deps "${J1}" "${J0B}" "${J0C}") "${DIR}/JOB_stage2B.sh" "${SAMPLES}");          echo "2B  annotate        : ${J2B}"
fi
if [ "${RUN[3]}" = 1 ]; then
    J3=$(sbatch --parsable $(deps "${J2B}" "${J2AR}") --array=1-"${N_FILTER_SETTINGS}" "${DIR}/JOB_stage3.sh" "${SAMPLES}")
    echo "3   filter (array)  : ${J3}"
fi

echo ""
echo "Submitted. Monitor with:  squeue -u \$USER"