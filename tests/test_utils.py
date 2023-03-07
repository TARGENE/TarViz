from tarviz.utils import load_pipeline_params, bQTL_list, http_variant_info, \
    http_ensemble_annotations
import os

NEXTFLOW_RUNDIR = os.path.join("tests", "nextflow_rundir")
PIPELINE_CONFIG_FILE = os.path.join(NEXTFLOW_RUNDIR, "nextflow.config")

def test_parse_pipeline_params():
    params = load_pipeline_params(PIPELINE_CONFIG_FILE)
    assert params == {
        "PARAMETER_PLAN": "FROM_ACTORS",
        "BQTLS": "data/bqtls.csv",
        "TRANS_ACTORS": "data/trans_acting_factors.csv",
        "EXTRA_COVARIATES": "data/extra_covariates.txt",
        "ORDERS": "1,2",
        "LD_BLOCKS": "data/VDR_LDblocks.txt",
        "NB_PCS": "6",
        "TRAITS_CONFIG": "data/ukbconfig.yaml",
        "ESTIMATORFILE": "data/estimator.yaml",
        "POSITIVITY_CONSTRAINT": "0.01",
        "PHENOTYPES_BATCH_SIZE": "0",
        "SIEVE_PVAL": "0.05",
        "NB_VAR_ESTIMATORS": "100",
        "MAX_TAU": "1.0",
        "OUTDIR": "$launchDir/results",
        "UKBB_BGEN_FILES": "/exports/igmm/eddie/UK-BioBank-53116/imputed/ukb_53116_chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}.{bgen,sample,bgen.bgi}",
        "UKBB_BED_FILES": "/exports/igmm/eddie/UK-BioBank-53116/genotypes/ukb_53116_chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}.{bed,bim,fam}",
        "QC_FILE": "/exports/igmm/eddie/UK-BioBank-53116/imputed/ukb_snp_qc.txt",
        "FLASHPCA_EXCLUSION_REGIONS": "data/exclusion_regions_hg19.txt",
        "WITHDRAWAL_LIST": "data/ukbb_participant_withdrawal.txt",
        "ENCRYPTED_DATASET": "/exports/igmm/eddie/UK-BioBank-53116/phenotypes/ukb45981.enc_ukb",
        "ENCODING_FILE": "/exports/igmm/eddie/UK-BioBank-53116/phenotypes/encoding.ukb"
    }

def test_bQTL_list():
    bqtls = bQTL_list(NEXTFLOW_RUNDIR)
    print(bqtls)
    assert bqtls == [
        "rs35405640", 
        "rs11544037", 
        "rs961320", 
        "rs12892296", 
        "rs10694113", 
        "rs11994337", 
        "rs13323956"
        ]
    
def test_http_variant_info():
    response = http_variant_info("rs35405640")
    assert response == {
        'phenotypes': [],
        'most_severe_consequence': 'intergenic_variant', 
        'mappings': [{
            'start': 150373493, 
            'allele_string': 'C/A/T', 
            'strand': 1, 
            'location': '3:150373493-150373493', 
            'end': 150373493, 
            'coord_system': 'chromosome', 
            'seq_region_name': '3', 
            'ancestral_allele': 'G', 
            'assembly_name': 'GRCh38'}], 
        'minor_allele': 'A', 
        'synonyms': [], 
        'name': 'rs35405640', 
        'MAF': 0.05351, 
        'evidence': ['Frequency', '1000Genomes', 'TOPMed', 'gnomAD'], 
        'ambiguity': 'H', 
        'var_class': 'SNP', 
        'source': 'Variants (including SNPs and indels) imported from dbSNP'
        }