DATA_COLUMNS = [
    "Column1","PARAMETER_TYPE","TREATMENTS","CASE","CONTROL",
    "TARGET","CONFOUNDERS","COVARIATES","INITIAL_ESTIMATE",
    "TMLE_ESTIMATE","TMLE_STD","TMLE_PVALUE","TMLE_LWB",
    "TMLE_UPB","ONESTEP_ESTIMATE","ONESTEP_STD","ONESTEP_PVALUE",
    "ONESTEP_LWB","ONESTEP_UPB","LOG","Treatment_1","Treatment_2",
    "Treatment_3","Case_1","Case_2","Case_3","Control_1","Control_2","Control_3",
    "WF_QVALUE"]

MT_ADJUSTEMENT_METHODS = ['TMLE_PVALUE', 'ONESTEP_PVALUE', "WF_QVALUE"]
#TRAIT_ADJUSTED_TMLE_PVALUE
ENSEMBL_URL = "https://rest.ensembl.org"

ANNOTATION_FEATURES = ["band", "gene", "transcript", "cds", "exon", "repeat", "simple", "misc", "variation", "somatic_variation", "structural_variation", "somatic_structural_variation", "constrained", "regulatory", "motif", "other_regulatory", "array_probe", "mane"]

PRIORITY_FEATURES = ["gene" ,"regulatory", "trait", "V2G"]