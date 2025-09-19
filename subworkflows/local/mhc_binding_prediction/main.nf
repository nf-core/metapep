nextflow.enable.dsl=2

// Bring in the modules this subworkflow orchestrates.
// NOTE: Keep these paths stable relative to projectDir.
include { PREPARE_PREDICTION_INPUT } from "${projectDir}/modules/local/prepare_prediction_input"
include { MHCFLURRY               }  from "${projectDir}/modules/local/mhcflurry"
include { MHCNUGGETS              }  from "${projectDir}/modules/local/mhcnuggets"



/* ──────────────────────────────────────────────────────────────────────────────
   POSTPROCESS for MHCflurry
   - Takes: (meta, input_csv) and (meta, predicted_csv)
   - Joins them safely (auto-detects swapped inputs), computes:
       * prediction_score  (IC50→binding affinity in [0,1])
       * rank_score        (presentation percentile / rank)
   - Emits a small, normalized TSV per input chunk for later merging.
   - Runs in a minimal pandas container (or conda env) only.
   ────────────────────────────────────────────────────────────────────────────── */
process MHCFLURRY_POSTPROCESS {
  label 'process_short'
  tag "${meta.id}"

  // Either conda OR container 
  conda "conda-forge::pandas=1.5.2"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
      'biocontainers/pandas:1.5.2' }"

  input:
  tuple val(meta), path(input_csv), path(predicted_csv)

  output:
  path "${meta.file_id}_mhcflurry_reduced.tsv", emit: reduced
  path "versions.yml",                            emit: versions

  shell:
  '''
  set -euo pipefail

  export FID="!{meta.file_id}"
  export INP="!{input_csv}"
  export PRED="!{predicted_csv}"

  python - <<'PY'
import os, sys, math
import pandas as pd

inp_path  = os.environ.get('INP',  "")
pred_path = os.environ.get('PRED', "")
fid       = os.environ.get('FID',  "output")

def read_csv(path):
    try:
        return pd.read_csv(path)
    except Exception as e:
        print(f"[ERROR] Could not read CSV: {path}", file=sys.stderr)
        print(str(e), file=sys.stderr)
        raise

inp  = read_csv(inp_path)
pred = read_csv(pred_path)

def looks_prepare(df):
    cols = {c.lower() for c in df.columns}
    return {"peptide","allele","peptide_id","allele_id"}.issubset(cols)

# Auto-fix swapped inputs: if needed, swap silently with a warning.
if not looks_prepare(inp) and looks_prepare(pred):
    print("[WARN] input_csv looks like predicted; swapping files.", file=sys.stderr)
    inp, pred = pred, inp

inp.columns  = [c.lower() for c in inp.columns]
pred.columns = [c.lower() for c in pred.columns]

# Input (prepared) must have peptide/allele plus IDs
need_inp = {"peptide","allele","peptide_id","allele_id"}
miss = need_inp - set(inp.columns)
if miss:
    raise KeyError(f"[ERROR] Missing columns in prepared CSV {inp_path}: {sorted(miss)}. Have: {list(inp.columns)}")

# MHCflurry outputs may vary by version; accept common column names.
affinity_col   = next((c for c in ("mhcflurry_affinity","affinity") if c in pred.columns), None)
percentile_col = next((c for c in ("mhcflurry_presentation_percentile","presentation_percentile","percentile","rank") if c in pred.columns), None)
if affinity_col is None:
    raise KeyError(f"[ERROR] No affinity column found in {pred_path}. Columns: {list(pred.columns)}")
if percentile_col is None:
    raise KeyError(f"[ERROR] No percentile/rank column found in {pred_path}. Columns: {list(pred.columns)}")

# Left-join by (peptide, allele). Validate: many predicted rows → one prepared row.
df = pred.merge(inp, on=["peptide","allele"], how="left", validate="many_to_one")

def ic50_to_ba(x):
    try:
        x = float(x)
    except Exception:
        return None
    if x > 50000:
        x = 50000.0
    # Map IC50 (nM) to binding affinity [0..1] with a log transform
    return 1.0 - (math.log10(x)/math.log10(50000.0))

# Unified scores
df["prediction_score"] = df[affinity_col].map(ic50_to_ba)
df["rank_score"]       = df[percentile_col]

# Coalesce IDs coming from either side of the merge (_x/_y)
def coalesce(df, *cols):
    s = None
    for c in cols:
        if c in df.columns:
            s = df[c] if s is None else s.combine_first(df[c])
    return s

pid = coalesce(df, "peptide_id_y", "peptide_id_x", "peptide_id")
aid = coalesce(df, "allele_id_y" , "allele_id_x" , "allele_id")


df["peptide_id_final"] = pid
df["allele_id_final"]  = aid

out_cols = ["peptide_id_final","allele_id_final","prediction_score","rank_score"]
missing_out = [c for c in out_cols if c not in df.columns]
if missing_out:
    raise KeyError("[ERROR] Missing columns after merge: "
                   f"{missing_out}. Available: {list(df.columns)}")

out = df[out_cols].rename(columns={"peptide_id_final":"peptide_id",
                                   "allele_id_final":"allele_id"})

out_path = fid + "_mhcflurry_reduced.tsv"
out.to_csv(out_path, sep="\t", index=False)

with open("versions.yml","w") as fh:
    fh.write("MHCFLURRY_POSTPROCESS:\\n")
    fh.write(f"  pandas: {pd.__version__}\\n")
PY
  '''
}


/* ──────────────────────────────────────────────────────────────────────────────
   POSTPROCESS for MHCnuggets
   - Similar role as MHCflurry postprocess; accepts slightly different column
     conventions (e.g., "hla" → "allele").
   ────────────────────────────────────────────────────────────────────────────── */
process MHCNUGGETS_POSTPROCESS {
  label 'process_short'
  tag "${meta.id}"

  conda "conda-forge::pandas=1.5.2"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
      'biocontainers/pandas:1.5.2' }"

  input:
  tuple val(meta), path(input_csv), path(predicted_csv)

  output:
  path "${meta.file_id}_mhcnuggets_reduced.tsv", emit: reduced
  path "versions.yml",                            emit: versions

  shell:
  '''
  set -euo pipefail
  export FID="!{meta.file_id}"
  export INP="!{input_csv}"
  export PRED="!{predicted_csv}"

  python - <<'PY'
import os, sys, math
import pandas as pd

inp  = pd.read_csv(os.environ['INP'])
pred = pd.read_csv(os.environ['PRED'])
fid  = os.environ.get('FID',"output")
STRICT = int(os.environ.get('STRICT_UNMAPPED','0'))

inp.columns  = [c.lower() for c in inp.columns]
pred.columns = [c.lower() for c in pred.columns]

need = {"peptide","allele","peptide_id","allele_id"}
if not need.issubset(set(inp.columns)):
    raise KeyError("Prepared CSV missing %s; have: %s" % (sorted(need - set(inp.columns)), list(inp.columns)))

# MHCnuggets sometimes uses 'hla' for allele; normalize to 'allele'
if "hla" in pred.columns and "allele" not in pred.columns:
    pred = pred.rename(columns={"hla":"allele"})

if not {"peptide","allele"}.issubset(set(pred.columns)):
    raise KeyError("Predictor output missing peptide/allele; have: %s" % list(pred.columns))

df = pred.merge(inp, on=["peptide","allele"], how="left", validate="many_to_one")

def ic50_to_ba(x):
    try: x=float(x)
    except: return None
    if x>50000: x=50000.0
    return 1.0 - (math.log10(x)/math.log10(50000.0))

affinity_col = next((c for c in ("ic50","affinity","nm","affinity_nm") if c in df.columns), None)
rank_col     = next((c for c in ("rank","percentile","rank_percentile") if c in df.columns), None)

df["prediction_score"] = df[affinity_col].map(ic50_to_ba) if affinity_col else None
df["rank_score"]       = df[rank_col] if rank_col else None

out = df[["peptide_id","allele_id","prediction_score","rank_score"]].copy()

# Write unmapped rows for inspection. Optionally fail in STRICT mode.
bad = out[out["peptide_id"].isna() | out["allele_id"].isna()]
if not bad.empty:
    bad.to_csv(f"{fid}_UNMAPPED_ROWS.tsv", sep="\t", index=False)
    if STRICT:
        raise SystemExit(f"Unmapped rows: {len(bad)} (strict mode)")
    out = out.dropna(subset=["peptide_id","allele_id"])

# Use compact integer dtypes where feasible
out["peptide_id"] = out["peptide_id"].astype("int64")
out["allele_id"]  = out["allele_id"].astype("int16")

out.to_csv(fid + "_mhcnuggets_reduced.tsv", sep="\t", index=False)

with open("versions.yml","w") as fh:
    fh.write("MHCNUGGETS_POSTPROCESS:\\n")
    fh.write("  pandas: %s\\n" % pd.__version__)
PY
  '''
}



/* ──────────────────────────────────────────────────────────────────────────────
   Merge many *reduced.tsv files into a single predictions.tsv.gz
   - Input: a single path value that is a LIST of TSV paths (collect() upstream)
   - the header is taken from the first file, then append the data rows from all
     subsequent files, and compress (bgzip if available, else gzip).
   ────────────────────────────────────────────────────────────────────────────── */
process MERGE_MHCFLURRY_REDUCED {
  label 'process_short'

  conda "conda-forge::pandas=1.5.2"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
      'biocontainers/pandas:1.5.2' }"

  input:
  // Expect exactly ONE emitted value that itself is the list of file paths
  path tsvs

  output:
  path "predictions.tsv.gz", emit: predictions
  path "versions.yml",      emit: versions

  script:
  """
  set -euo pipefail

  # Nextflow → Bash: turn "file1 file2 file3" into a Bash array
  FILES=( ${tsvs.join(' ')} )

  if [ \${#FILES[@]} -eq 0 ]; then
    echo "[ERROR] No TSVs provided to merge." >&2
    exit 1
  fi

  # Header from the first, then data rows from each file
  if command -v bgzip >/dev/null 2>&1; then
    { head -n1 "\${FILES[0]}"; for f in "\${FILES[@]}"; do tail -n +2 "\$f"; done; } | bgzip -c > predictions.tsv.gz
  else
    { head -n1 "\${FILES[0]}"; for f in "\${FILES[@]}"; do tail -n +2 "\$f"; done; } | gzip  -c > predictions.tsv.gz
  fi

  {
    echo "MERGE_MHCFLURRY_REDUCED:"
    echo "  python: \$(python --version | sed 's/Python //g')"
    echo "  pandas: \$(python -c 'import pandas as p; print(p.__version__)')"
  } > versions.yml
  """
}


/* ──────────────────────────────────────────────────────────────────────────────
   Subworkflow: MHC_BINDING_PREDICTION
   - Input: channel of peptide chunks, each as (meta, TSV/CSV)
     * meta must carry at least: meta.id and (ideally) meta.file_id
   - Steps:
     1) PREPARE_PREDICTION_INPUT normalizes & expands inputs for the chosen tool
     2) Run the requested predictor (mhcflurry | mhcnuggets )
     3) POSTPROCESS converts predictor outputs into normalized reduced TSVs
     4) MERGE concatenates all reduced tables into predictions.tsv.gz
   - Outputs:
     * predictions : single gzipped TSV with normalized columns
     * versions    : combined versions from all steps
   ────────────────────────────────────────────────────────────────────────────── */
workflow MHC_BINDING_PREDICTION {

  take:
    ch_in_peptides   // Tuples: (meta, tsv) — chunks incl. meta.alleles etc.

  main:
    // every tuple has a stable file_id for naming outputs
    def ch_peptides_safe = ch_in_peptides.map { m, f ->
      def fid = m.file_id ?: "${m.id ?: f.baseName}_${f.baseName}"
      [ m + [file_id: fid], f ]
    }

    // Prepare inputs for the chosen tool; emits (meta, prepared_csv)
    def prep = PREPARE_PREDICTION_INPUT(ch_peptides_safe)

    def tool = (params.pred_method ?: '').toLowerCase().trim()

    def out_predictions = Channel.empty()
    def out_versions    = prep.versions

    if (tool == 'mhcflurry') {
      // MHCFLURRY module should emit (meta, predicted_csv)
      def pred   = MHCFLURRY( prep.prepared )

      // Align (meta, predicted_csv) with (meta, prepared_csv) by meta
      def joined = pred.predicted.join(prep.prepared).map { m, pred_csv, in_csv -> tuple(m, in_csv, pred_csv) }

      // Postprocess to normalized reduced TSVs, then merge
      def post   = MHCFLURRY_POSTPROCESS( joined )
      def merged = MERGE_MHCFLURRY_REDUCED( post.reduced.collect() )

      out_predictions = merged.predictions
      out_versions    = out_versions.mix(pred.versions).mix(post.versions).mix(merged.versions)
    }
    else if (tool == 'mhcnuggets-class-1') {
      // If your MHCNUGGETS module expects (meta, tsv), we pass the prepared CSV as "tsv".
      def pred   = MHCNUGGETS( prep.prepared.map { m, csv -> tuple(m, csv) } )
      def joined = pred.predicted.join(prep.prepared).map { m, pred_csv, in_csv -> tuple(m, in_csv, pred_csv) }
      def post   = MHCNUGGETS_POSTPROCESS( joined )
      def merged = MERGE_MHCFLURRY_REDUCED( post.reduced.collect() )

      out_predictions = merged.predictions
      out_versions    = out_versions.mix(pred.versions).mix(post.versions).mix(merged.versions)
    }
    else {
      error "Unsupported params.pred_method='${params.pred_method}'. Allowed: mhcflurry | mhcnuggets-class-1 "
    }

  emit:
    predictions = out_predictions
    versions    = out_versions
}
