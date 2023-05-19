DATA_COLUMNS = [
    'PARAMETER_TYPE', 'TREATMENTS', 'TARGET', 'TRAIT_ADJUSTED_TMLE_PVALUE',
    'TMLE_PVALUE', 'TMLE_ESTIMATE', 'TMLE_LWB', 'TMLE_UPB', 'TMLE_STD', 
    'ONESTEP_PVALUE', 'ONESTEP_ESTIMATE', 'ONESTEP_LWB', 'ONESTEP_UPB', 'ONESTEP_STD',
    'INITIAL_ESTIMATE', 'CASE', 'CONTROL', 'CONFOUNDERS', 'COVARIATES', 'LOG']

MT_ADJUSTEMENT_METHODS = ['TRAIT_ADJUSTED_TMLE_PVALUE', 'TMLE_PVALUE', 'ONESTEP_PVALUE']

ENSEMBL_URL = "https://rest.ensembl.org"

ANNOTATION_FEATURES = ["band", "gene", "transcript", "cds", "exon", "repeat", "simple", "misc", "variation", "somatic_variation", "structural_variation", "somatic_structural_variation", "constrained", "regulatory", "motif", "other_regulatory", "array_probe", "mane"]