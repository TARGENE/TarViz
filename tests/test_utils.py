from tarviz.utils import load_pipeline_params, bQTLs_data, http_variant_info, \
    http_ensemble_annotations, http_ensembl_binding_matrix
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

def test_bQTLs_data():
    bqtls = bQTLs_data(NEXTFLOW_RUNDIR)
    assert bqtls.ID.tolist() == [
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
    assert 'phenotypes' in response
    assert 'mappings' in response
    assert response["name"] == "rs35405640"
    
def test_http_ensemble_annotations():
    response = http_ensemble_annotations('3', 150373493, 150373493)
    assert len(response) > 1

def test_http_ensembl_binding_matrix():
    response = http_ensembl_binding_matrix("ENSPFM0487")
    assert "elements" in response
    assert "stable_id" in response

