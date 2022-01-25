process PREDICT_EPITOPES {
    label 'process_low'
    label 'cache_lenient'

    // TODO: conda
    conda (params.enable_conda ? { exit 1 "Conda is currently not available for metapep" } : null)
    container 'skrakau/metapep:dev'

    input:
    path(peptides)

    output:
    path "*predictions.tsv",            emit:   ch_epitope_predictions
    path "*prediction_warnings.log",    emit:   ch_epitope_prediction_warnings
    path "versions.yml",                emit:   versions

    script:
    def pred_method           = params.pred_method
    switch (pred_method) {
        case "syfpeithi":
            pred_method_version = "1.0";
            break;
        case "mhcflurry":
            pred_method_version = "1.4.3";
            break;
        case "mhcnuggets-class-1":
            pred_method_version = "2.3.2";
            break;
        case "mhcnuggets-class-2":
            pred_method_version = "2.3.2";
            break;
        }
    """

    # Extract allele name from file header
    allele_name="\$(head -n1 "$peptides" | fgrep '#' | cut -f2 -d'#')"
    allele_id="\$(head -n1 "$peptides" | fgrep '#' | cut -f3 -d'#')"

    out_basename="\$(basename "$peptides" .txt)"
    out_predictions="\$out_basename"_predictions.tsv
    out_warnings="\$out_basename"_prediction_warnings.log

    # Create output header
    echo "peptide_id	prediction_score	allele_id" >"\$out_predictions"

    # Process file
    # The --syfpeithi-norm flag enables score normalization when syfpeithi is
    # used and is ignored otherwise
    if ! epytope_predict.py --peptides "$peptides" \\
                    --method "$pred_method" \\
                    --method_version "$pred_method_version" \\
                    --syfpeithi-norm \\
                    "\$allele_name" \\
                    2>stderr.log \\
                    | tail -n +2 \\
                    | cut -f 1,3 \\
                    | sed -e "s/\$/	\$allele_id/" \\
                    >>"\$out_basename"_predictions.tsv; then
        cat stderr.log >&2
        exit 1
    fi

    # Filter stderr for warnings and pass them on in the warnings channel
    fgrep WARNING stderr.log  | sort -u >"\$out_warnings" || :

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        fred2: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('Fred2').version)")
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
        pyvcf: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pyvcf').version)")
        mhcflurry: \$(mhcflurry-predict --version 2>&1 | sed 's/^mhcflurry //; s/ .*\$//')
        mhcnuggets: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('mhcnuggets').version)")
    END_VERSIONS
    """
}
