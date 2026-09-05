#!/bin/bash
# Run nhmmer for each key in a multi-model HMM file, convert tblouts to BED,
# and merge to a single BED for use with `krep mask --hmm-bed`.
#
# Usage:
#   run_nhmmer_layer.sh <hmm_file> <target_hmmer_db> <out_prefix> [keys...]
#
# Environment:
#   E       E-value threshold (default 1e-2)
#   CPU     threads per nhmmer job (default 2)
#
# Example:
#   run_nhmmer_layer.sh Dfam_curatedonly.hmm chm13.hmmerdb hmm_chr1 \
#                       MIR L2 CR1_Mam

set -e
HMM_FILE=$1
TARGET_DB=$2
OUT_PREFIX=$3
shift 3
KEYS=${*:-"MIR L2 CR1_Mam"}
E=${E:-1e-2}
CPU=${CPU:-2}

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

: > "${OUT_PREFIX}.log"
for K in $KEYS; do
    echo "=== $K ===" >> "${OUT_PREFIX}.log"
    hmmfetch "$HMM_FILE" "$K" > "${OUT_PREFIX}_${K}.hmm" 2>/dev/null
    nhmmer --noali --tblout "${OUT_PREFIX}_${K}.tbl" -E "$E" --cpu "$CPU" \
           "${OUT_PREFIX}_${K}.hmm" "$TARGET_DB" >/dev/null 2>&1
    python3 "${SCRIPT_DIR}/nhmmer_to_bed.py" "${OUT_PREFIX}_${K}.tbl" "${OUT_PREFIX}_${K}.bed"
    echo "$K hits: $(wc -l < "${OUT_PREFIX}_${K}.bed")" >> "${OUT_PREFIX}.log"
done

cat "${OUT_PREFIX}"_*.bed | sort -k1,1 -k2,2n | awk 'BEGIN{OFS="\t"}
    {if (NR==1 || $1!=c || $2>e){if(NR>1)print c,s,e,f; c=$1;s=$2;e=$3;f=$4} else {e=$3}}
    END{print c,s,e,f}' > "${OUT_PREFIX}_merged.bed"

echo "merged: $(wc -l < "${OUT_PREFIX}_merged.bed")" >> "${OUT_PREFIX}.log"
echo "Wrote ${OUT_PREFIX}_merged.bed"
