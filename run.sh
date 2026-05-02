#!/usr/bin/env bash
# Full pipeline: STRIDE -> tripeptide extraction -> angle computation -> plot
# Usage: bash run.sh [pdb_directory] [target_residue]

set -euo pipefail

PDB_DIR="${1:-data/pdb_files}"
RESIDUE="${2:-ARG}"
CHUNK=1000

export PDB_PATH="$PDB_DIR"

mkdir -p tmp/pdb tmp/stride tmp/windows output

# enumerate structures
echo "scanning $PDB_DIR ..."
find "$PDB_DIR" -name "*.pdb.gz" -exec basename {} .pdb.gz \; > tmp/all_ids.txt
TOTAL=$(wc -l < tmp/all_ids.txt | tr -d ' ')
echo "$TOTAL structures found"

# --- phase 1: STRIDE ---
echo "phase 1: secondary structure assignment"
chunk_idx=0
while IFS= read -r id; do
    echo "$id" >> "tmp/chunk_${chunk_idx}.txt"
    n=$(wc -l < "tmp/chunk_${chunk_idx}.txt" | tr -d ' ')
    if [ "$n" -ge "$CHUNK" ]; then
        snakemake -s pipeline/stride.smk --cores 4 --keep-going \
            --config batch_file="tmp/chunk_${chunk_idx}.txt" pdb_path="$PDB_DIR" 2>/dev/null
        rm "tmp/chunk_${chunk_idx}.txt"
        chunk_idx=$((chunk_idx + 1))
    fi
done < tmp/all_ids.txt

# process remaining
if [ -f "tmp/chunk_${chunk_idx}.txt" ]; then
    snakemake -s pipeline/stride.smk --cores 4 --keep-going \
        --config batch_file="tmp/chunk_${chunk_idx}.txt" pdb_path="$PDB_DIR" 2>/dev/null
    rm "tmp/chunk_${chunk_idx}.txt"
fi

DONE=$(find tmp/stride/ -name "*.stride" 2>/dev/null | wc -l | tr -d ' ')
echo "  completed: $DONE / $TOTAL"

# --- phase 2: extract tripeptide windows ---
echo "phase 2: tripeptide window extraction"
processed=0
for ss_file in tmp/stride/*.stride; do
    [ -f "$ss_file" ] || continue
    pdb_id=$(basename "$ss_file" .stride)
    dst="tmp/windows/${RESIDUE}_${pdb_id}.tsv"
    [ -f "$dst" ] && continue
    python3 src/stride_parser.py "$ss_file" "$dst" "$RESIDUE" 2>/dev/null || true
    processed=$((processed + 1))
    [ $((processed % 5000)) -eq 0 ] && echo "  $processed extracted..."
done
echo "  windows: $(find tmp/windows/ -name '*.tsv' | wc -l | tr -d ' ')"

# --- phase 3: angle computation ---
echo "phase 3: vector angle computation"
python3 src/compute_vectors.py tmp/windows "$RESIDUE" output/angles.tsv

# --- phase 4: visualization ---
echo "phase 4: generating plot"
python3 src/visualize.py

rm -f tmp/all_ids.txt
echo "pipeline complete!"
