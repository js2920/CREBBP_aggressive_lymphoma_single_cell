#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Batch Cell Ranger GEX alignment for the Shabanas 10x dataset.
# - Discovers all SLX-*.SIG*.tar archives under ${data_root}
# - Extracts FASTQs per barcode (once, unless force_reextract=true)
# - Maps barcodes to friendly sample names via sample_attribution.txt
# - Runs `cellranger count` sequentially, reserving 28 cores per job
# -----------------------------------------------------------------------------

set -euo pipefail
shopt -s nullglob

############################ USER SETTINGS ####################################
data_root="/home/gusti/CREBBP_aggressive_lymphoma_single_cell/data/raw_fastq"
sample_attr="${data_root}/sample_attribution.txt"
fastq_workspace="${data_root}/cellranger_fastqs"
run_root="${data_root}/cellranger_runs"

ref_dir="/home/gusti/CREBBP_aggressive_lymphoma_single_cell/data/refdata-gex-GRCm39-2024-A"
cellranger_bin="cellranger"

cores=22
mem_gb=128
create_bam=false          # true keeps BAMs, false saves space
force_reextract=false     # set true to re-untar even if FASTQs already exist
###############################################################################

mkdir -p "$fastq_workspace" "$run_root"

command -v python3 >/dev/null 2>&1 || { echo "❌ python3 is required."; exit 1; }
[[ -x "$cellranger_bin" ]] || { echo "❌ Cell Ranger binary not found → $cellranger_bin"; exit 1; }
[[ -d "$ref_dir"       ]] || { echo "❌ Reference folder not found → $ref_dir"; exit 1; }
[[ -d "$data_root"     ]] || { echo "❌ FASTQ archive root missing → $data_root"; exit 1; }
[[ -f "$sample_attr"   ]] || { echo "⚠️ Sample attribution sheet missing → $sample_attr"; }

clean_id() {
  local raw="${1:-sample}"
  raw="${raw// /_}"
  raw=$(printf '%s' "$raw" | sed -E 's/[^A-Za-z0-9_]+/_/g' | sed -E 's/_+/_/g' | sed -E 's/^_+//; s/_+$//')
  [[ -z "$raw" ]] && raw="sample"
  printf '%s' "$raw"
}

metadata_lines=$(python3 - "$sample_attr" <<'PY'
import sys, re, pathlib
path = pathlib.Path(sys.argv[1])
if not path.exists():
    sys.exit(0)
emit = False
for raw in path.read_text().splitlines():
    if not raw.strip():
        continue
    if raw.startswith("For metadata"):
        emit = True
        continue
    if not emit or raw.startswith("ID\t"):
        continue
    parts = raw.split('\t')
    if len(parts) < 5:
        continue
    barcode = parts[-1].strip()
    sample = parts[-2].strip()
    clean = re.sub(r'[^A-Za-z0-9_]+', '_', sample.replace(' ', '_')).strip('_')
    if not clean:
        clean = barcode
    print(f"{barcode}\t{clean}")
PY
)

declare -A barcode_to_label=()
if [[ -n "${metadata_lines//[[:space:]]/}" ]]; then
  while IFS=$'\t' read -r barcode label; do
    [[ -z "$barcode" ]] && continue
    barcode_to_label["$barcode"]="$label"
  done <<< "$metadata_lines"
else
  echo "⚠️ No metadata rows parsed from $sample_attr; will label runs with barcodes."
fi

declare -A barcode_tar_map=()
while IFS= read -r -d '' tarball; do
  fname=$(basename "$tarball")
  barcode=$(awk -F'.' '{print $2}' <<<"$fname")
  [[ -z "$barcode" ]] && continue
  barcode_tar_map["$barcode"]+="$tarball"$'\n'
done < <(find "$data_root" -type f -name "SLX-*.SIG*.tar" ! -name "*lostreads*" -print0)

if [[ ${#barcode_tar_map[@]} -eq 0 ]]; then
  echo "❌ No SLX-*.SIG*.tar archives found under $data_root"
  exit 1
fi

echo "🔎 Found ${#barcode_tar_map[@]} barcode group(s) to process."
mapfile -t sorted_barcodes < <(printf "%s\n" "${!barcode_tar_map[@]}" | sort)

pushd "$run_root" >/dev/null

for barcode in "${sorted_barcodes[@]}"; do
  tar_entries="${barcode_tar_map[$barcode]}"
  mapfile -t tar_paths < <(printf '%s' "$tar_entries" | sed '/^$/d')
  [[ ${#tar_paths[@]} -eq 0 ]] && continue

  sample_label="${barcode_to_label[$barcode]:-$barcode}"
  safe_label=$(clean_id "$sample_label")
  run_id="${barcode}_${safe_label}_GEX"
  fastq_dest="${fastq_workspace}/${barcode}"

  echo "================================================================"
  echo "📂 Barcode: ${barcode}"
  echo "   Label  : ${sample_label}"
  echo "   FASTQs :"
  printf '     • %s\n' "${tar_paths[@]}"

  needs_extract=true
  if [[ "$force_reextract" == false && -d "$fastq_dest" ]]; then
    existing_r1=("${fastq_dest}"/*_R1_*.fastq.gz)
    if [[ ${#existing_r1[@]} -gt 0 ]]; then
      needs_extract=false
    fi
  fi

  if [[ "$needs_extract" == true ]]; then
    rm -rf "$fastq_dest"
    mkdir -p "$fastq_dest"
    for tar_path in "${tar_paths[@]}"; do
      echo "   ↪ Extracting $(basename "$tar_path")"
      tar -xf "$tar_path" -C "$fastq_dest"
    done
    printf '%s\n' "${tar_paths[@]}" > "${fastq_dest}/.source_tarballs.txt"
  else
    echo "   ↪ FASTQs already extracted in $fastq_dest (set force_reextract=true to refresh)"
  fi

  r1_files=("${fastq_dest}"/*_R1_*.fastq.gz)
  if [[ ${#r1_files[@]} -eq 0 ]]; then
    echo "❌ No R1 FASTQs detected for ${barcode}; skipping."
    continue
  fi

  if [[ -d "${run_root}/${run_id}/outs" ]]; then
    echo "✅ Existing Cell Ranger output detected → ${run_root}/${run_id}/outs (skipping)."
    continue
  fi

  echo "▶ Running Cell Ranger: --id=${run_id}, --sample=${barcode}"
  "$cellranger_bin" count \
    --id="$run_id" \
    --fastqs="$fastq_dest" \
    --sample="$barcode" \
    --transcriptome="$ref_dir" \
    --create-bam="$create_bam" \
    --nosecondary \
    --localcores="$cores" \
    --localmem="$mem_gb"

  status=$?
  if [[ $status -eq 0 ]]; then
    echo "✅ ${run_id} finished → ${run_root}/${run_id}/outs"
  else
    echo "❌ ${run_id} failed with exit code ${status}"
  fi
done

popd >/dev/null
echo "🎉 Alignment loop finished."



