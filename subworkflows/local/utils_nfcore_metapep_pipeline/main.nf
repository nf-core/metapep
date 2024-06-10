//
// Subworkflow with functionality specific to the nf-core/metapep pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFVALIDATION_PLUGIN } from '../../nf-core/utils_nfvalidation_plugin'
include { paramsSummaryMap          } from 'plugin/nf-validation'
include { fromSamplesheet           } from 'plugin/nf-validation'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { dashedLine                } from '../../nf-core/utils_nfcore_pipeline'
include { nfCoreLogo                } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { workflowCitation          } from '../../nf-core/utils_nfcore_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    pre_help_text = nfCoreLogo(monochrome_logs)
    post_help_text = '\n' + workflowCitation() + '\n' + dashedLine(monochrome_logs)
    def String workflow_command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"
    UTILS_NFVALIDATION_PLUGIN (
        help,
        workflow_command,
        pre_help_text,
        post_help_text,
        validate_params,
        "nextflow_schema.json"
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )
    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //
    ch_samplesheet = params.show_supported_models ? Channel.empty() : Channel.fromPath(input, checkIfExists: true)

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs, multiqc_report.toList())
        }

        completionSummary(monochrome_logs)

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    // Exit if peptide length parameters are exchanged
    if (params.min_pep_len > params.max_pep_len) {
        error "The minimum peptide length needs to be smaller or equal than the maximum. See 'https://nf-co.re/metapep/dev/parameters' for more information."
    }
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def citation_text = [
            "Tools used in the workflow included:",
            "Entrez (Maglott et al. 2005),",
            "Prodigal (Hyatt et al. 2010),",
            "Python (Python Core Team 2022),",
            "R (R Core Team 2022),",
            "Pandas (The pandas development team 2022),",
            "Epytope (FRED2) (Schubert et al. 2016),",
            params["pred_method"]=="syfpeithi" ? "SYFPEITHI (Rammensee et al. 1999)" : "",
            params["pred_method"]=="mhcflurry" ? "MHCFlurry (O'Donnell et al. 2020)" : "",
            params["pred_method"]=="mhcnuggets-class-1" ? "Shao (Rammensee et al. 2020)" : "",
            params["pred_method"]=="mhcnuggets-class-2" ? "Shao (Rammensee et al. 2020)" : "",
            "and",
            "MultiQC (Ewels et al. 2016).",
            "For a more detailed list check the references and tool versions",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    def reference_text = [
            "<li>Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., & Notredame, C. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316-319. doi: <a href='https://doi.org/10.1038/nbt.3820'>10.1038/nbt.3820</a></li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: <a href='https://doi.org/10.1093/bioinformatics/btw354'>/10.1093/bioinformatics/btw354</li>",
            "<li>Ewels, P. A., Peltzer, A., Fillinger, S., Patel, H., Alneberg, J., Wilm, A., Garcia, M. U., Di Tommaso, P., & Nahnsen, S. (2020). The nf-core framework for community-curated bioinformatics pipelines. Nature Biotechnology, 38(3), 276-278. doi: <a href='https://doi.org/10.1038/s41587-020-0439-x'>10.1038/s41587-020-0439-x</a></li>",
            "<li>Grüning, B., Dale, R., Sjödin, A., Chapman, B. A., Rowe, J., Tomkins-Tinch, C. H., Valieris, R., Köster, J., & Bioconda Team. (2018). Bioconda: sustainable and comprehensive software distribution for the life sciences. Nature Methods, 15(7), 475–476. doi: <a href='https://doi.org/10.1038/s41592-018-0046-7'>10.1038/s41592-018-0046-7</a></li>",
            "<li>Hyatt, D., Chen, GL., LoCascio, P. F., Land M. L., Larimer F. W., Hauser L. J. (2010). Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics 11, 119. doi: <a href='https://doi.org/10.1186/1471-2105-11-119'>10.1186/1471-2105-11-119</a></li>",
            "<li>Maglott D, Ostell J, Pruitt KD, Tatusova T. (2005) Entrez Gene: gene-centered information at NCBI. Nucleic Acids Res. 2005 Jan 1;33(Database issue):D54-8. Update in: Nucleic Acids Res. 2007 Jan;35(Database issue):D26-31. doi: <a href='https://doi.org/10.1093/nar/gki031'>10.1093/nar/gki031</a>.</li>",
            "<li>O'Donnell T. J., Rubinsteyn A., Laserson U., (2020). MHCflurry 2.0: Improved Pan-Allele Prediction of MHC Class I-Presented Peptides by Incorporating Antigen Processing. Cell Systems 11, 42-48. doi: <a href='https://doi.org/10.1016/j.cels.2020.06.010'>10.1016/j.cels.2020.06.010</a>.</li>",
            "<li>Python Core Team (2022). Python: A dynamic, open source programming language. Python Software Foundation. <a href='https://www.python.org/'>https://www.python.org/</a>.</li>",
            "<li>Rammensee H., Bachmann J., Emmerich N. P., Bachor O. A., Stevanović S. (1999). SYFPEITHI: database for MHC ligands and peptide motifs. Immunogenetics 1999 Nov;50(3-4):213-9. doi: <a href='https://doi.org/10.1007/s002510050595'>10.1007/s002510050595</a>.</li>",
            "<li>R Core Team (2022). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. <a href='https://www.R-project.org/'>https://www.R-project.org/</a>.</li>",
            "<li>Schubert, B., Walzer, M., Brachvogel, H-P., Sozolek, A., Mohr, C., and Kohlbacher, O. (2016). FRED 2 - An Immunoinformatics Framework for Python. Bioinformatics 2016. doi: <a href='https://doi.org/10.1093/bioinformatics/btw113'>10.1093/bioinformatics/btw113</a>.</li>",
            "<li>Shao X. M., Bhattacharya R., Huang J., Sivakumar I. K. A., Tokheim C., Zheng L., Hirsch D., Kaminow B., Omdahl A., Bonsack M., Riemer A. B., Velculescu V. E., Anagnostou V., Pagel K. A., Karchin R. (2020). High-Throughput Prediction of MHC Class I and II Neoantigens with MHCnuggets. Cancer Immunol Res. 2020 Mar;8(3):396-408. doi: <a href='https://doi.org/10.1158/2326-6066.cir-19-0464'>10.1158/2326-6066.cir-19-0464</a>.</li>",
            "<li>The pandas development team. (2022). pandas-dev/pandas: Pandas (v1.5.2). Zenodo. doi: <a href='https://doi.org/10.5281/zenodo.7344967'>10.5281/zenodo.7344967</a>.</li>",
            "<li>da Veiga Leprevost, F., Grüning, B. A., Alves Aflitos, S., Röst, H. L., Uszkoreit, J., Barsnes, H., Vaudel, M., Moreno, P., Gatto, L., Weber, J., Bai, M., Jimenez, R. C., Sachsenberg, T., Pfeuffer, J., Vera Alvarez, R., Griss, J., Nesvizhskii, A. I., & Perez-Riverol, Y. (2017). BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics (Oxford, England), 33(16), 2580–2582. doi: <a href='https://doi.org/10.1093/bioinformatics/btx192'>10.1093/bioinformatics/btx192</a></li>"
    ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        String[] manifest_doi = meta.manifest_map.doi.tokenize(",")
        for (String doi_ref: manifest_doi) temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
