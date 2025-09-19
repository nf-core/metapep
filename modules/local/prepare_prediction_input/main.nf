process PREPARE_PREDICTION_INPUT {
  label 'process_short'
  tag "${meta.id}"
  container "https://depot.galaxyproject.org/singularity/mhcgnomes:1.8.6--pyh7cba7a3_0"

  input:
  tuple val(meta), path(tsv)

  output:
  tuple val(meta), path("${meta.file_id}_mhcflurry_input.csv"), emit: prepared
  path  "versions.yml",                                           emit: versions

  shell:
  '''
  set -euo pipefail
  
  FID="!{meta.file_id}"
  TSV="!{tsv}"

  # --- 1) Header line: "#<allele_name>#<allele_id>"
  header_line=$(head -n1 "$TSV")
  if ! printf '%s\n' "$header_line" | grep -q '^#'; then
    echo "[ERROR] First line does not contain a comment header with allele (#...#...)" >&2
    exit 1
  fi
  allele_name=$(printf '%s' "$header_line" | awk -F'#' '{print $2}')
  allele_id=$(printf  '%s' "$header_line" | awk -F'#' '{print $3}')
  if [ -z "$allele_name" ] || [ -z "$allele_id" ]; then
    echo "[ERROR] Could not read allele_name/allele_id from '#name#id': $header_line" >&2
    exit 1
  fi

  # --- 2) Drop comment lines, normalize delimiters, unify header names
  CLEAN_TSV="${FID}_clean.tsv"
  awk '
    BEGIN{ first=1 }
    /^#/ { next }  # skip comment lines
    {
      line=$0
      sub(/^[ \t]+/,"",line); sub(/[ \t]+$/,"",line)
      gsub(/[ \t]+/,"\t",line)
      $0=line
    }
    first==1 {
      n=split($0, H, "\t")
      for(i=1;i<=n;i++){
        h=H[i]; sub(/^[ \t]+/,"",h); sub(/[ \t]+$/,"",h)
        hl=tolower(h)
        if(hl=="peptide_id" || hl=="id" || hl=="peptideid" || hl=="peptide-id") H[i]="peptide_id"
        else if(hl=="peptide_sequence" || hl=="sequence" || hl=="peptide" || hl=="seq") H[i]="peptide_sequence"
      }
      out=H[1]; for(i=2;i<=n;i++) out=out"\t"H[i]
      print out
      first=0; next
    }
    { print }
  ' "$TSV" > "$CLEAN_TSV"

  # --- 3) Validate header (peptide_id + peptide_sequence)
  hdr=$(head -n1 "$CLEAN_TSV")
  if ! printf '%s\n' "$hdr" | grep -Eq '^peptide_id\t'; then
    echo "[ERROR] Expected column peptide_id in ${CLEAN_TSV}. Header was:" >&2
    echo "  $hdr" >&2
    exit 1
  fi
  if ! printf '%s\n' "$hdr" | grep -Eq '\tpeptide_sequence(\t|$)'; then
    echo "[ERROR] Expected column peptide_sequence in ${CLEAN_TSV}. Header was:" >&2
    echo "  $hdr" >&2
    exit 1
  fi

  # --- 4) CSV for MHCflurry: peptide,allele,peptide_id,allele_id
  outcsv="${FID}_mhcflurry_input.csv"
  {
    echo "peptide,allele,peptide_id,allele_id"
    awk -F'\t' -v OFS=',' -v a="$allele_name" -v aid="$allele_id" '
      NR==1 {
        for(i=1;i<=NF;i++){ if($i=="peptide_sequence") p=i; if($i=="peptide_id") id=i }
        if(!p || !id){ print "[ERROR] Header without peptide_sequence and/or peptide_id. Header was:", $0 > "/dev/stderr"; exit 1 }
        next
      }
      { if(NF>0) print $p, a, $id, aid }
    ' "$CLEAN_TSV"
  } > "$outcsv"

  # --- 5) Ensure CSV contains no TABs
  if grep -q $'\t' "$outcsv"; then
    echo "[ERROR] Output contains TABs; should be pure CSV: $outcsv" >&2
    head -n3 "$outcsv" | cat -A >&2
    exit 1
  fi

  # --- 6) Minimal versions.yml
  {
    echo "NFCORE_METAPEP:METAPEP:MHC_BINDING_PREDICTION:PREPARE_PREDICTION_INPUT:"
    echo "  python: NA"
    echo "  pandas: NA"
    echo "  mhcgnomes: NA"
  } > versions.yml
  '''
}
