#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/metapep
========================================================================================
 nf-core/metapep Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/metapep
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    // TODO nf-core: Add to this help message with new command line parameters
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/metapep --input '*_R{1,2}.fastq.gz' -profile docker

    Mandatory arguments:
      --input [file]                  Path to input TSV file containing: condition, type, path, alleles. Must contain a corresponding header.
      -profile [str]                  Configuration profile to use. Can use multiple (comma separated)
                                      Available: conda, docker, singularity, test, awsbatch, <institute> and more

    Options:
      --genome [str]                  Name of iGenomes reference
      --prodigal_mode [str]           Prodigal mode, 'meta' or 'single'. Default: 'meta'.
      --ncbi_key [str]                NCBI key for faster download from Entrez databases.
      --ncbi_email [str]              Email address for NCBI Entrez database access. Required if downloading proteins from NCBI.
      --min_pep_len [int]             Min. peptide length to generate.
      --max_pep_len [int]             Max. peptide length to generate.
      --sample_n [int]                Number of peptides to subsample for each condition. Default: false
      --pred_method [str]             Epitope prediction method to use. One of [syfpeithi, mhcflurry, mhcnuggets-class-1, mhcnuggets-class-2]. Default: syfpeithi.
      --pred_chunk_size               Maximum chunk size (#peptides) for epitope prediction jobs

    Other options:
      --outdir [file]                 The output directory where the results will be saved
      --publish_dir_mode [str]        Mode for publishing results in the output directory. Available: symlink, rellink, link, copy, copyNoFollow, move (Default: copy)
      --email [email]                 Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --email_on_fail [email]         Same as --email, except only send mail if the workflow is not successful
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    AWSBatch options:
      --awsqueue [str]                The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion [str]               The AWS Region for your AWS Batch job to run on
      --awscli [str]                  Path to the AWS CLI tool
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)

/*
 * Check and process input file
 */
if (!params.input)
    exit 1, "Missing input. Please specify --input."
if (!hasExtension(params.input, "tsv"))
    exit 1, "Input file specified with --input must have a '.tsv' extension."

switch (params.pred_method) {
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
    default:
        exit 1, "Epitope prediction method specified with --pred_method not recognized."
}

// TODO for 'proteins' type
// - allow weight input for type 'proteins' as well! (for now use equal weight ?)

// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Input']            = params.input
summary['Prodigal mode']    = params.prodigal_mode
summary['Min. peptide length']   = params.min_pep_len
summary['Max. peptide length']   = params.max_pep_len
summary['Peptide Subsampling'] = params.sample_n ? "$params.sample_n per condition" : "disabled"
summary['Prediction method']     = params.pred_method
summary['Prediction chunk size'] = params.pred_chunk_size
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-metapep-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/metapep Workflow Summary'
    section_href: 'https://github.com/nf-core/metapep'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf(".csv") > 0) filename
                      else null
                }

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file "software_versions.csv"

    script:
    // TODO nf-core: Get all tools to print their version number here
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/*
 * Create database tables for input
 */
process create_db_tables {
    publishDir "${params.outdir}/db_tables", mode: params.publish_dir_mode,
        saveAs: {filename -> "$filename" }
    input:
    file input_file from Channel.value(file(params.input))

    output:
    file "microbiomes.tsv" into ch_microbiomes                    // microbiome_id, microbiome_path, microbiome_type, weights_path
    file "conditions.tsv" into ch_conditions                      // condition_id, condition_name, microbiome_id
    file "alleles.tsv"  into ch_alleles                           // allele_id, allele_name
    file "conditions_alleles.tsv" into ch_conditions_alleles      // condition_id, allele_id

    script:
    """
    create_db_tables.py -i ${input_file} \
                        -m microbiomes.tsv \
                        -c conditions.tsv \
                        -a alleles.tsv \
                        -ca conditions_alleles.tsv
    """
}

// ####################################################################################################

ch_microbiomes
    // Read microbiomes table
    .splitCsv(sep:'\t', header:true)
    // Convert paths to files
    .map {
        row ->
        row.microbiome_path = file(row.microbiome_path, checkIfExists: true)
        row
    }
    // Split into types
    .branch {
        row->
        taxa:      row.microbiome_type == 'taxa'
        proteins : row.microbiome_type == 'proteins'
        assembly:  row.microbiome_type == 'assembly'
        bins:      row.microbiome_type == 'bins'
        other:     true
    }
    .set{ch_microbiomes_branch}

// Emit warning about unknown data types
ch_microbiomes_branch.other.map { row -> log.warn("Ignoring row in input sheet: Unknown type '${row.microbiome_type}'") }

// TAXA
ch_microbiomes_branch.taxa
    .multiMap { row ->
            ids: row.microbiome_id
            files: row.microbiome_path
        }
    .set { ch_taxa_input }

// PROTEINS
ch_microbiomes_branch.proteins
    .multiMap { row ->
            ids: row.microbiome_id
            files: row.microbiome_path
        }
    .set { ch_proteins_input }

// ASSEMBLY
ch_microbiomes_branch.assembly
    .multiMap { row ->
            ids: row.microbiome_id
            files: row.microbiome_path
            bin_basenames: false
        }
    .set { ch_assembly_input }

// BINS
ch_microbiomes_branch.bins
    .branch {
            row ->
            folders : row.microbiome_path.isDirectory()
            archives : row.microbiome_path.isFile()
            other: true
        }
    .set{ ch_microbiomes_bins }

// The file ending we expect for FASTA files
fasta_suffix = ~/(?i)[.]fa(sta)?(.gz)?$/

// BINS - LOCAL FOLDERS
ch_microbiomes_bins.folders
    .multiMap { row ->
        def bin_files = row.microbiome_path.listFiles().findAll{ it.name =~ fasta_suffix }
        ids           : Collections.nCopies((int) bin_files.size(), row.microbiome_id)
        files         : bin_files
        bin_basenames : bin_files.collect{ it.name - fasta_suffix }
    }.set { ch_microbiomes_bins_folders }

// BINS - LOCAL OR REMOTE ARCHIVES
ch_microbiomes_bins.archives
    .multiMap { row ->
        ids : row.microbiome_id
        files: row.microbiome_path
    }
    .set{ch_microbiomes_archives}

process unpack_bin_archives {
    input:
    val microbiome_id from ch_microbiomes_archives.ids
    path microbiome_path from ch_microbiomes_archives.files

    output:
    tuple val(microbiome_id), file('unpacked/*') into ch_microbiomes_unpacked_archives

    script:
    """
    mkdir -v unpacked
    tar -C unpacked -vxf "$microbiome_path"
    """
}

ch_microbiomes_bins_archives = Channel.empty()
ch_microbiomes_unpacked_archives
    .multiMap { microbiome_id, bin_files ->
        bin_files = bin_files.findAll{ it.name =~ fasta_suffix }
        ids           : Collections.nCopies((int) bin_files.size(), microbiome_id)
        files         : bin_files
        bin_basenames : bin_files.collect{ it.name - fasta_suffix }
    }.set{ch_microbiomes_bins_archives}

// Concatenate the channels for nucleotide based inputs
ch_nucl_input_ids           = ch_assembly_input.ids.concat(ch_microbiomes_bins_archives.ids.flatten(), ch_microbiomes_bins_folders.ids.flatten())
ch_nucl_input_files         = ch_assembly_input.files.concat(ch_microbiomes_bins_archives.files.flatten(), ch_microbiomes_bins_folders.files.flatten())
ch_nucl_input_bin_basenames = ch_assembly_input.bin_basenames.concat(ch_microbiomes_bins_archives.bin_basenames.flatten(), ch_microbiomes_bins_folders.bin_basenames.flatten())


// ####################################################################################################

ch_weights = Channel.empty()
ch_microbiomes
    .splitCsv(sep:'\t', skip: 1)
    .map { microbiome_id, microbiome_path, microbiome_type, weights_path ->
            if (microbiome_type != 'taxa' && weights_path) [microbiome_id, weights_path]
        }
    .view().multiMap { microbiome_id, weights_path ->
            microbiome_ids: microbiome_id
            weights_paths: weights_path
        }.set { ch_weights }

/*
 * Download proteins from entrez
 */
process download_proteins {
    publishDir "${params.outdir}", mode: params.publish_dir_mode,
        saveAs: {filename ->
                    if (filename.indexOf(".fasta.gz") == -1) "entrez_data/$filename"
                    else null
        }

    input:
    val    microbiome_ids     from   ch_taxa_input.ids.collect()
    file   microbiome_files   from   ch_taxa_input.files.collect()

    output:
    file   "proteins.entrez.tsv.gz"            into   ch_entrez_proteins
    file   "taxa_assemblies.tsv"               into   ch_entrez_assemblies
    file   "entities_proteins.entrez.tsv"      into   ch_entrez_entities_proteins  // protein_tmp_id (accessionVersion), entity_name (taxon_id)
    file   "microbiomes_entities.entrez.tsv"   into   ch_entrez_microbiomes_entities  // entity_name, microbiome_id, entity_weight

    script:
    def key = params.ncbi_key
    def email = params.ncbi_email
    def microbiome_ids = microbiome_ids.join(' ')
    """
    # provide new home dir to avoid permission errors with Docker and other artefacts
    export HOME="\${PWD}/HOME"
    download_proteins_entrez.py --email $email \
                                --key $key \
                                -t $microbiome_files \
                                -m $microbiome_ids \
                                -p proteins.entrez.tsv.gz \
                                -ta taxa_assemblies.tsv \
                                -ep entities_proteins.entrez.tsv \
                                -me microbiomes_entities.entrez.tsv
    """
}

/*
 * Predict proteins from contigs
 */
process predict_proteins {
    publishDir "${params.outdir}/prodigal", mode: params.publish_dir_mode,
        saveAs: {filename ->
                    if (filename.indexOf(".fasta") == -1) "$filename"
                    else null
        }

    input:
    val microbiome_id from ch_nucl_input_ids
    val bin_basename from ch_nucl_input_bin_basenames
    file microbiome_file from ch_nucl_input_files

    output:
    val microbiome_id into ch_pred_proteins_microbiome_ids                  // Emit microbiome ID
    val bin_basename into ch_pred_proteins_bin_basename
    file("proteins.pred_${microbiome_id}*.tsv.gz") into ch_pred_proteins     // Emit protein tsv
    file "coords.pred_${microbiome_id}*.gff"

    script:
    def mode = params.prodigal_mode
    def name = bin_basename ? "${microbiome_id}.${bin_basename}" : "${microbiome_id}"
    """
    gzip -c -d $microbiome_file | prodigal \
                -f gff \
                -o coords.pred_${name}.gff \
                -a proteins.pred_${name}.fasta \
                -p $mode

    echo -e "protein_tmp_id\tprotein_sequence" > proteins.pred_${name}.tsv
    fasta_to_tsv.py --remove-asterisk --input proteins.pred_${name}.fasta >> proteins.pred_${name}.tsv
    gzip proteins.pred_${name}.tsv
    """
}

/*
 * Assign entity weights for input type 'assembly' and 'bins'
 */
process assign_nucl_entity_weights {
    publishDir "${params.outdir}/db_tables", mode: params.publish_dir_mode,
        saveAs: {filename -> "$filename" }

    input:
    val  microbiome_ids     from  ch_weights.microbiome_ids.collect().ifEmpty([])
    path weights_files      from  ch_weights.weights_paths.collect().ifEmpty([])

    output:
    path   "microbiomes_entities.nucl.tsv"    into   ch_nucl_microbiomes_entities  // entity_name, microbiome_id, entity_weight

    script:
    microbiome_ids = microbiome_ids.join(' ')
    """
    assign_entity_weights.py \
        --microbiome-ids $microbiome_ids \
        --weights-files $weights_files \
        --out microbiomes_entities.nucl.tsv
    """
}

/*
 * concat files and assign new, unique ids for all proteins (from different sources)
 */
process generate_protein_and_entity_ids {
    publishDir "${params.outdir}/db_tables", mode: params.publish_dir_mode,
        saveAs: {filename -> "$filename" }

    input:
    // Predicted Proteins
    path   predicted_proteins                  from       ch_pred_proteins.collect().ifEmpty([])
    val    predicted_proteins_microbiome_ids   from       ch_pred_proteins_microbiome_ids.collect().ifEmpty([])
    val    predicted_proteins_bin_basenames    from       ch_pred_proteins_bin_basename.collect().ifEmpty([])
    // Entrez Proteins
    path   entrez_proteins                     from       ch_entrez_proteins.ifEmpty([])
    path   entrez_entities_proteins            from       ch_entrez_entities_proteins.ifEmpty([])       //   protein_tmp_id (accessionVersion), entity_name (taxon_id)
    path   entrez_microbiomes_entities         from       ch_entrez_microbiomes_entities.ifEmpty([])    //   entity_name, microbiome_id, entity_weight
    // Bare Proteins
    path   bare_proteins                       from       ch_proteins_input.files.collect().ifEmpty([])
    path   bare_proteins_microbiome_ids        from       ch_proteins_input.ids.collect().ifEmpty([])

    output:
    path   "proteins.tsv.gz"                        into   ch_proteins
    path   "entities_proteins.tsv"                  into   ch_entities_proteins
    path   "entities.tsv"                           into   ch_entities
    path   "microbiomes_entities.no_weights.tsv"    into   ch_microbiomes_entities_noweights  // microbiome_id, entitiy_id  (no weights yet!)

    script:
    predicted_proteins_microbiome_ids = predicted_proteins_microbiome_ids.join(' ')
    predicted_proteins_bin_basenames  = predicted_proteins_bin_basenames.collect{ it ? it : "__ISASSEMBLY__" }.join(' ')
    """
    generate_protein_and_entity_ids.py \
        --predicted-proteins                  $predicted_proteins                  \
        --predicted-proteins-microbiome-ids   $predicted_proteins_microbiome_ids   \
        --predicted-proteins-bin-basenames    $predicted_proteins_bin_basenames    \
        --entrez-proteins                     "$entrez_proteins"                   \
        --entrez-entities-proteins            "$entrez_entities_proteins"          \
        --entrez-microbiomes-entities         "$entrez_microbiomes_entities"       \
        --bare-proteins                       $bare_proteins                       \
        --bare-proteins-microbiome-ids        $bare_proteins_microbiome_ids        \
        --out-proteins                        proteins.tsv.gz                      \
        --out-entities-proteins               entities_proteins.tsv                \
        --out-entities                        entities.tsv                         \
        --out-microbiomes-entities            microbiomes_entities.no_weights.tsv
    """
}

/*
 * Create microbiome_entities
 */
process finalize_microbiome_entities {
    publishDir "${params.outdir}/db_tables", mode: params.publish_dir_mode,
        saveAs: {filename -> "$filename" }

    input:
    path   entrez_microbiomes_entities        from       ch_entrez_microbiomes_entities.ifEmpty([])
    path   nucl_microbiomes_entities          from       ch_nucl_microbiomes_entities.ifEmpty([])
    path   microbiomes_entities_noweights     from       ch_microbiomes_entities_noweights
    path   entities                           from       ch_entities

    output:
    path   "microbiomes_entities.tsv"    into   ch_microbiomes_entities  // entity_id, microbiome_id, entity_weight

    script:

    """
    finalize_microbiome_entities.py \
        -eme $entrez_microbiomes_entities \
        -nme $nucl_microbiomes_entities \
        -menw $microbiomes_entities_noweights \
        -ent "$entities" \
        -o microbiomes_entities.tsv
    """
    // TODO add checking if for microbiome_id either no weight or weights for all entities are given
}

/*
 * Generate peptides
 */
process generate_peptides {
    publishDir "${params.outdir}", mode: params.publish_dir_mode,
        saveAs: {filename -> "db_tables/$filename" }

    input:
    file proteins from ch_proteins

    output:
    file "peptides.tsv.gz" into ch_peptides                // peptide_id, peptide_sequence
    file "proteins_peptides.tsv" into ch_proteins_peptides // protein_id, peptide_id, count
    //file "proteins_lengths.tsv"

    script:
    def min_pep_len = params.min_pep_len
    def max_pep_len = params.max_pep_len
    """
    generate_peptides.py -i $proteins \
                         -min $min_pep_len \
                         -max $max_pep_len \
                         -p "peptides.tsv.gz" \
                         -pp "proteins_peptides.tsv" \
                         -l "proteins_lengths.tsv"
    """
}

/*
 * Split prediction tasks (peptide, allele) into chunks of peptides that are to
 * be predicted against the same allele for parallel prediction
 */
process split_pred_tasks {
    input:
    path  peptides              from  ch_peptides
    path  proteins_peptides     from  ch_proteins_peptides
    path  entities_proteins     from  ch_entities_proteins
    path  microbiomes_entities  from  ch_microbiomes_entities
    path  conditions            from  ch_conditions
    path  conditions_alleles    from  ch_conditions_alleles
    path  alleles               from  ch_alleles
    // The tables are joined to map peptide -> protein -> microbiome -> condition -> allele
    // and thus to enumerate, which (peptide, allele) combinations have to be predicted.

    output:
    path "peptides_*.txt" into ch_epitope_prediction_chunks

    script:
    def pred_chunk_size       = params.pred_chunk_size
    def subsampling = params.sample_n ? "--sample_n ${params.sample_n}" : ""
    """
    gen_prediction_chunks.py --peptides "$peptides" \
                             --protein-peptide-occ "$proteins_peptides" \
                             --entities-proteins-occ "$entities_proteins" \
                             --microbiomes-entities-occ "$microbiomes_entities" \
                             --conditions "$conditions" \
                             --condition-allele-map "$conditions_alleles" \
                             --max-chunk-size $pred_chunk_size \
                             $subsampling \
                             --alleles "$alleles" \
                             --outdir .
    """
}

/*
 * Perform epitope prediction
 */
process predict_epitopes {
    input:
    path peptides from ch_epitope_prediction_chunks.flatten()

    output:
    path "*predictions.tsv" into ch_epitope_predictions
    path "*prediction_warnings.log" into ch_epitope_prediction_warnings

    script:
    def pred_method           = params.pred_method
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
    if ! epytope_predict.py --peptides "$peptides" \
                       --method "$pred_method" \
                       --method_version "$pred_method_version" \
		       --syfpeithi-norm \
                       "\$allele_name" \
                       2>stderr.log \
                       | tail -n +2 \
                       | cut -f 1,3 \
                       | sed -e "s/\$/	\$allele_id/" \
                       >>"\$out_basename"_predictions.tsv; then
        cat stderr.log >&2
        exit 1
    fi

    # Filter stderr for warnings and pass them on in the warnings channel
    fgrep WARNING stderr.log  | sort -u >"\$out_warnings" || :
    """
}

/*
 * Merge prediction results from peptide chunks into one prediction result
 */
 // gather chunks of predictions and merge them already to avoid too many input files for `merge_predictions` process
 // (causing "sbatch: error: Batch job submission failed: Pathname of a file, directory or other parameter too long")
process merge_predictions_buffer {

    input:
    path predictions from ch_epitope_predictions.buffer(size: 1000, remainder: true)
    path prediction_warnings from ch_epitope_prediction_warnings.buffer(size: 1000, remainder: true)

    output:
    path "predictions.buffer_*.tsv" into ch_predictions_buffer
    path "prediction_warnings.buffer_*.log" into ch_prediction_warnings_buffer

    script:
    def single = predictions instanceof Path ? 1 : predictions.size()
    def merge = (single == 1) ? 'cat' : 'csvtk concat -t'
    """
    [[ ${predictions[0]} =~  peptides_(.*)_predictions.tsv ]];
    uname="\${BASH_REMATCH[1]}"
    echo \$uname

    $merge $predictions > predictions.buffer_\$uname.tsv
    sort -u $prediction_warnings > prediction_warnings.buffer_\$uname.log
    """
}

process merge_predictions {
    publishDir "${params.outdir}", mode: params.publish_dir_mode,
        saveAs: {filename -> filename.endsWith(".log") ? "logs/$filename" : "db_tables/$filename"}

    input:
    path predictions from ch_predictions_buffer.collect()
    path prediction_warnings from ch_prediction_warnings_buffer.collect()

    output:
    path "predictions.tsv.gz" into ch_predictions
    path "prediction_warnings.log"

    script:
    def single = predictions instanceof Path ? 1 : predictions.size()
    def merge = (single == 1) ? 'cat' : 'csvtk concat -t'
    """
    $merge $predictions | gzip > predictions.tsv.gz
    sort -u $prediction_warnings > prediction_warnings.log
    """
}

/*
 * Generate figures
 */

process prepare_score_distribution {
    publishDir "${params.outdir}/figures/prediction_scores", mode: params.publish_dir_mode

    input:
    file predictions from ch_predictions
    file proteins_peptides from ch_proteins_peptides
    file entities_proteins from ch_entities_proteins
    file microbiomes_entities from ch_microbiomes_entities
    file conditions from  ch_conditions
    file conditions_alleles from  ch_conditions_alleles
    file alleles from ch_alleles

    output:
    file "prediction_scores.allele_*.tsv" into ch_prep_prediction_scores

    script:
    """
    prepare_score_distribution.py --predictions "$predictions" \
                            --protein-peptide-occ "$proteins_peptides" \
                            --entities-proteins-occ "$entities_proteins" \
                            --microbiomes-entities-occ "$microbiomes_entities" \
                            --conditions "$conditions" \
                            --condition-allele-map "$conditions_alleles" \
                            --alleles "$alleles" \
                            --outdir .
    """
}

process plot_score_distribution {
    publishDir "${params.outdir}/figures", mode: params.publish_dir_mode

    input:
    file prep_scores from ch_prep_prediction_scores.flatten()
    file alleles from ch_alleles
    file conditions from ch_conditions

    output:
    file "prediction_score_distribution.*.pdf"

    script:
    """
    [[ ${prep_scores} =~ prediction_scores.allele_(.*).tsv ]];
    allele_id="\${BASH_REMATCH[1]}"
    echo \$allele_id

    plot_score_distribution.R --scores $prep_scores \
                                   --alleles $alleles \
                                   --conditions $conditions \
                                   --allele_id \$allele_id \
                                   --method ${params.pred_method}
    """
}


process prepare_entity_binding_ratios {
    publishDir "${params.outdir}/figures/entity_binding_ratios", mode: params.publish_dir_mode

    input:
    file predictions from ch_predictions
    file proteins_peptides from ch_proteins_peptides
    file entities_proteins from ch_entities_proteins
    file microbiomes_entities from ch_microbiomes_entities
    file conditions from  ch_conditions
    file conditions_alleles from  ch_conditions_alleles
    file alleles from ch_alleles

    output:
    file "entity_binding_ratios.allele_*.tsv" into ch_prep_entity_binding_ratios

    script:
    """
    prepare_entity_binding_ratios.py --predictions "$predictions" \
                            --protein-peptide-occ "$proteins_peptides" \
                            --entities-proteins-occ "$entities_proteins" \
                            --microbiomes-entities-occ "$microbiomes_entities" \
                            --conditions "$conditions" \
                            --condition-allele-map "$conditions_alleles" \
                            --alleles "$alleles" \
                            --method ${params.pred_method} \
                            --outdir .
    """
}

process plot_entity_binding_ratios {
    publishDir "${params.outdir}/figures", mode: params.publish_dir_mode

    input:
    file prep_entity_binding_ratios from ch_prep_entity_binding_ratios.flatten()
    file alleles from ch_alleles

    output:
    file "entity_binding_ratios.*.pdf"

    script:
    """
    [[ ${prep_entity_binding_ratios} =~ entity_binding_ratios.allele_(.*).tsv ]];
    allele_id="\${BASH_REMATCH[1]}"
    echo \$allele_id

    plot_entity_binding_ratios.R --binding-rates $prep_entity_binding_ratios \
                                   --alleles $alleles \
                                   --allele_id \$allele_id
    """
}


/*
 * Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    file output_docs from ch_output_docs
    file images from ch_output_docs_images

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/metapep] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/metapep] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    def mqc_report = null

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$projectDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$projectDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$projectDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$projectDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/metapep] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            mail_cmd.execute() << email_html
            log.info "[nf-core/metapep] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/metapep]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/metapep]${c_red} Pipeline completed with errors${c_reset}-"
    }

}

def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/metapep v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}
