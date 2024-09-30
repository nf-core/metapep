process PREDICT_EPITOPES {
    label 'process_low'
    label 'process_long'
    label 'error_retry'

    conda "bioconda::epytope=3.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/epytope:3.3.1--pyh7cba7a3_0' :
        'biocontainers/epytope:3.3.1--pyh7cba7a3_0' }"

    input:
    path(peptides)

    output:
    path "*predictions.tsv",            emit:   ch_epitope_predictions
    path "*pred_warnings.log",          emit:   ch_epitope_prediction_warnings
    path "versions.yml",                emit:   versions

    script:
    """
    # create folder for MHCflurry downloads to avoid permission problems when running pipeline with docker profile and mhcflurry selected
    mkdir -p mhcflurry-data
    export MHCFLURRY_DATA_DIR=./mhcflurry-data
    # specify MHCflurry release for which to download models, need to be updated here as well when MHCflurry will be updated
    export MHCFLURRY_DOWNLOADS_CURRENT_RELEASE=1.4.0

    # Extract software versions from container
    mhcflurry_version=\$(mhcflurry-predict --version 2>&1 | sed 's/^mhcflurry //; s/ .*\$//')
    mhcnuggets_version=\$(python -c "import pkg_resources; print(pkg_resources.get_distribution('mhcnuggets').version)")

    # Syfpeithi is not an external software, but rather a matrix on which scoring is based on -> titled version 1.0 in epytope
    syfpeithi_version=1.0

    # Assign version based on method
    case $params.pred_method in

        "syfpeithi")
        pred_method_version=\$syfpeithi_version
        ;;

        "mhcflurry")
        pred_method_version=\$mhcflurry_version
        ;;

        "mhcnuggets-class-1" | "mhcnuggets-class-2")
        pred_method_version=\$mhcnuggets_version
        ;;
    esac

    # Extract allele name from file header
    allele_name="\$(head -n1 "$peptides" | fgrep '#' | cut -f2 -d'#')"
    allele_id="\$(head -n1 "$peptides" | fgrep '#' | cut -f3 -d'#')"

    out_basename="\$(basename "$peptides" .txt)"
    out_predictions="\$out_basename"_predictions.tsv
    out_warnings="\$out_basename"_pred_warnings.log

    # Create output header
    echo "peptide_id	prediction_score	allele_id" >"\$out_predictions"

    # Process file
    # The --syfpeithi-norm flag enables score normalization when syfpeithi is
    # used and is ignored otherwise
    if ! epytope_predict.py --peptides "$peptides" \\
                    --method "$params.pred_method" \\
                    --method_version "\$pred_method_version" \\
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
        epytope: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('epytope').version)")
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
        mhcflurry: \$mhcflurry_version
        mhcnuggets: \$mhcnuggets_version
        syfpeithi: \$syfpeithi_version
    END_VERSIONS
    """
}
