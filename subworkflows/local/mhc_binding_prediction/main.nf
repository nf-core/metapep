nextflow.enable.dsl=2

include { PREPARE_PREDICTION_INPUT } from "${projectDir}/modules/local/prepare_prediction_input"
include { MHCFLURRY               }  from "${projectDir}/modules/local/mhcflurry"
include { MHCNUGGETS              }  from "${projectDir}/modules/local/mhcnuggets"
include { POSTPROCESS_PREDICTIONS }  from "${projectDir}/modules/local/postprocess_prediction"
include { MERGE_POSTPROCESS       }  from "${projectDir}/modules/local/merge_postprocess"

workflow MHC_BINDING_PREDICTION {
  take:
    ch_in_peptides
    ch_alleles_file

  main:
  def SUPJSON_VAL = Channel.value( file("$projectDir/assets/supported_alleles.json") ) 
  def ALLELES_VAL = ch_alleles_file.map { file(it) }.first() 

  //prepare predictor-specific inputs (and id map for mhcnuggets)
  PREPARE_PREDICTION_INPUT( ch_in_peptides, SUPJSON_VAL, ALLELES_VAL )
  
  def ch_flurry_in  = PREPARE_PREDICTION_INPUT.out.flurry .map { m, flurry_csv, allele_txt -> tuple(m, flurry_csv) }
  def ch_nuggets_in = PREPARE_PREDICTION_INPUT.out.nuggets.map { m, nuggets_tsv, allele_txt ->
    def als = file(allele_txt).text.trim()
    def cls = (params.pred_method == 'mhcnuggets-class-2') ? 'II' : 'I'
    tuple(m + [ alleles_supported: als, mhc_class: cls ], nuggets_tsv)
  }

  def ch_idmap = PREPARE_PREDICTION_INPUT.out.idmap //for mhcnuggets, columns: peptide, allele, peptide_id, allele_id
  def out_vers = PREPARE_PREDICTION_INPUT.out.versions

  
  def ch_post_in 

  if (params.pred_method == 'mhcflurry') {
    MHCFLURRY(ch_flurry_in)
    
    //pair MHCflurry predictions with their input for post-processing
    ch_post_in = MHCFLURRY.out.predicted
      .map  { m, pred_csv  -> tuple(m.id, [m, pred_csv]) }
      .join ( ch_flurry_in.map { m, input_csv -> tuple(m.id, [m, input_csv]) } )
      .map  { id, L, R ->
        def (m1, pred_csv ) = L
        def (m2, input_csv) = R
        assert m1.id == m2.id
        // method-Label für Postprocess
        tuple(m1, 'mhcflurry', input_csv, pred_csv)
      }

    out_vers = out_vers.mix(MHCFLURRY.out.versions)

  } else if (params.pred_method in ['mhcnuggets-class-1','mhcnuggets-class-2']) {
    MHCNUGGETS( ch_nuggets_in )

    //pair MHCnuggets predictions with ID map for post-processing (to recover peptide/allele IDs)
    ch_post_in = MHCNUGGETS.out.predicted
      .map  { m, pred_csv  -> tuple(m.id, [m, pred_csv]) }
      .join ( ch_idmap.map  { m, idmap_csv -> tuple(m.id, [m, idmap_csv]) } )
      .map  { id, L, R ->
        def (m1, pred_csv ) = L
        def (m2, idmap_csv) = R
        assert m1.id == m2.id

        tuple(m1, 'mhcnuggets', idmap_csv, pred_csv)
      }

    out_vers = out_vers.mix(MHCNUGGETS.out.versions)

  } else {
    error "Unsupported --pred_method='${params.pred_method}'. Supported: mhcflurry | mhcnuggets-class-1 | mhcnuggets-class-2"
  }


  POSTPROCESS_PREDICTIONS( ch_post_in )
  out_vers = out_vers.mix(POSTPROCESS_PREDICTIONS.out.versions)

 //merge all reduced TSVs into predcitions.tsv.gz
  MERGE_POSTPROCESS( POSTPROCESS_PREDICTIONS.out.reduced.collect() )
  out_vers = out_vers.mix(MERGE_POSTPROCESS.out.versions)

  emit:
    predictions = MERGE_POSTPROCESS.out.predictions
    versions    = out_vers
}

