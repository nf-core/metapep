#!/usr/bin/env python3
import argparse, math, os, sys
import pandas as pd

def ic50_to_binding_affinity(x):
    """Convert IC50 value to binding affinity score between 0 and 1."""
    if x > 50000:
        x = 50000.0
    return 1.0 - (math.log10(x) / math.log10(50000.0))



def normalize_cols(df):
    """Convert all column names to lowercase."""
    df.columns = [c.lower() for c in df.columns]
    return df


def postprocess(method, inp_path, pred_path, out_path, strict_unmapped=False):

    ''' in the following peptide = peptide sequence and allele = allele name !

        mhcflurry inp columns (mhcflurry_input.csv): peptide,allele,peptide_id,allele_id 
        mhcflurry pred columns (mhcflurry_predicted.csv): peptide,allele,mhcflurry_affinity,mhcflurry_presentation_score,mhcflurry_presentation_percentile
        mhcnuggets inp columns (mhcnuggets_idmap.csv): peptide,allele,peptide_id,allele_id
        mhcnuggets pred columns (mhcnuggets_predicted.csv): peptide,ic50,rank,allele
        '''
    inp  = pd.read_csv(inp_path)
    pred = pd.read_csv(pred_path)
    inp  = normalize_cols(inp) 
    pred = normalize_cols(pred) 

   
    # Merge prediction results with input data based on peptide sequence and allele name
    # This links predictions back to original peptide/allele IDs
    df = pred.merge(inp, on=["peptide","allele"], how="left", suffixes=('_pred', '_inp'), validate="many_to_one")

    if method == "mhcflurry":
        # Mhcflurry uses mhcflurry_affinity (IC50-like) and presentation percentile
       aff_col = "mhcflurry_affinity" if "mhcflurry_affinity" in df.columns else None
       rank_col = "mhcflurry_presentation_percentile" if "mhcflurry_presentation_percentile" in df.columns else None
    
    else:  
       # Mhcnuggets uses ic50 values and human proteome rank
       aff_col = "ic50" if "ic50" in df.columns else None
       rank_col = "human_proteome_rank" if "human_proteome_rank" in df.columns else None


    df["prediction_score"] = df[aff_col].map(ic50_to_binding_affinity)
    df["rank_score"]       = df[rank_col] if rank_col is not None else pd.NA
    

   # Extract peptide and allele IDs from the original input data (not predictions)
   # mhcnuggets uses peptide_id/allele_id, mhcflurry uses peptide_id_inp/allele_id_inp
    peptide_id = "peptide_id_inp" if "peptide_id_inp" in df.columns else "peptide_id" 
    allele_id = "allele_id_inp" if "allele_id_inp" in df.columns else "allele_id"  

    # Create final output DataFrame with standardized column structure
    # Final format: peptide_id, allele_id, prediction_score, rank_score
    out = (df.rename(columns={peptide_id: "peptide_id", allele_id: "allele_id"})
        [["peptide_id","allele_id","prediction_score","rank_score"]])
    
    out.to_csv(out_path, sep="\t", index=False)
    print(f"Successfully created: {out_path} :)", file=sys.stderr)
    
   
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--method", choices=["mhcflurry","mhcnuggets"], required=True)
    ap.add_argument("--input",  required=True) #mhcflurry/mhcnuggets_input.csv
    ap.add_argument("--pred",   required=True) #mhcflurry/mhcnuggets_predicted.csv
    ap.add_argument("--out",    required=True) #mhcflurry/mhcnuggets_reduced.tsv

    args = ap.parse_args() 
    postprocess(args.method, args.input, args.pred, args.out)


if __name__ == "__main__":
    sys.exit(main())