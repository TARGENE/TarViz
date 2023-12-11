import pandas as pd
import os
import requests
import tables
import streamlit as st
from tarviz.constants import ENSEMBL_URL, DATA_COLUMNS



def result_file(nextflow_rundir):
    return os.path.join(nextflow_rundir, "preliminary_results.csv")

def raw_data_file(nextflow_rundir):
    return os.path.join(nextflow_rundir, "results", "tmle_inputs", "final.data.csv")

def load_pipeline_params(config_file):
    in_param_section = False
    params = {}
    for line in open(config_file, 'r').readlines():
        sline = line.strip()
        if sline and in_param_section:
            # End of the param section
            if sline == "}":
                return params
            # Skip Comments in the param file
            elif sline.startswith("//"):
                continue
            # Split key values and append dict
            else:
                key, val = sline.split("=")
                params[key.strip()] = val.strip().strip("'").strip('"')
        elif sline.startswith("params"):
            in_param_section = True

def get_bql(treatment_str, bqtl_set):
    for t in treatment_str.split("_&_"):
        if t in bqtl_set:
            return t
    return None

@st.cache_data
def load_data(bqtl_data=None, results_file=None):
    data = pd.read_csv(os.path.join(st.session_state["nextflow_rundir"], "results", results_file))
    if bqtl_data is not None:
         bqtl_set = set(bqtl_data.ID.unique())
         data["BQTL"] = [get_bql(x, bqtl_set) for x in data.TREATMENTS]
    return data 

@st.cache_data
def load_ref_alt_counts(ref_alt_data):
     ref_alt_data = pd.read_csv(os.path.join(st.session_state["nextflow_rundir"], "data", ref_alt_data))[DATA_COLUMNS]
     return ref_alt_data
    

@st.cache_data
def bQTLs_data(config_filename="nextflow.config"):
    params = load_pipeline_params(os.path.join(st.session_state.nextflow_rundir, config_filename))

    bqtls_filepath = os.path.join(st.session_state.nextflow_rundir, params["BQTLS"])

    if os.path.exists(bqtls_filepath):
        return pd.read_csv(bqtls_filepath, sep='\s+')
    return None

def http_variant_info(rsid):
    url = "".join((
        ENSEMBL_URL,
        "/variation/human/",
        rsid,
        "?",
        "phenotypes=1"
    ))
    return requests.get(url, headers={ "Content-Type": "application/json"}).json()

def http_ensemble_annotations(
        chr, start, end, 
        distance=100, 
        features=("gene", "regulatory", "motif")
        ):
    url = "".join((
        ENSEMBL_URL,
        "/overlap/region/human/", 
        chr, 
        ":", 
        str(start - distance), 
        "-", 
        str(end + distance),
        "?", 
        ";".join("".join(("feature=", f)) for f in features)
    ))
    return requests.get(url, headers={"Content-Type": "application/json"}).json()

def http_ensembl_binding_matrix(stable_matrix_id):
    url = "".join((
        ENSEMBL_URL,
        "/species/homo_sapiens/binding_matrix/", 
        stable_matrix_id, 
        "?", 
        "unit=frequencies", 
    ))
    return requests.get(url, headers={"Content-Type": "application/json"}).json()

@st.cache_data
def bqtls_hit_counts(data):
    counts = data.groupby(["BQTL"])["TARGET"].nunique().reset_index(name="COUNTS")
    counts.sort_values("COUNTS", ascending=False, inplace=True)

    return counts["BQTL"] + " (" + counts["COUNTS"].astype(str) + " hits)"

@st.cache_data
def feature_columns(df, feature):
    if feature == "motif":
        return df[["transcription_factor_complex", "binding_matrix_stable_id", "score", "start", "end", "strand"]]
    else:
        return df


### EXTRA ANNOTATION/EXPRESSION INFO 
@st.cache_data

def look_up_variant_gtex_tissue(variant_code_b38_ref_alt):
    
    string_variant = "chr"+variant_code_b38_ref_alt+"_b38"
    gtex = pd.read_hdf(st.session_state["gtex_file"])

    gtex[["VARIANT_CODE","ENSEMBL"]] = gtex["RSID"].str.split(',', expand=True)
    gtex_og= gtex[gtex["VARIANT_CODE"]== string_variant]
    st.dataframe(gtex_og[["ENSEMBL","#STUDY","PVALUE_FE", "BETA_FE", "STD_FE", "PVALUE_RE", "BETA_RE", "STD_RE",
            "PVALUE_RE2","STAT1_RE2", "STAT2_RE2", "PVALUE_BE", "I_SQUARE", "Q", "PVALUE_Q", "TAU_SQUARE"]])

    gtex_og= gtex[gtex["VARIANT_CODE"]== string_variant]

    gtex.columns = pd.MultiIndex.from_frame(
        pd.DataFrame(
            gtex
            .columns
            .str
            .split("_", n=1)
            .to_list()
        )
    )

    # accessing tissue
    tissue =  gtex.columns.get_level_values(1)
    return gtex_og, tissue

def look_up_variant_gtex(gtex, variant_code_b38_ref_alt, selection):
    gtex.columns = pd.MultiIndex.from_frame(
        pd.DataFrame(
            gtex
            .columns
            .str
            .split("_", n=1)
            .to_list()
        )
    )

    # accessing tissue
    tissue =  gtex.columns.get_level_values(1)

    return gtex.loc[:, tissue == selection]

@st.cache_data
# Open Targets V2G
def open_targets_v2g(variantId):
    """Function to look up vairant in targets and output V2G score in JSON format."""
    query_string = """
    query v2g($variantId: String!) {
      genesForVariant(variantId: $variantId) {
        gene {
          id
        }
        qtls {
        typeId
        sourceId
        aggregatedScore
        }
        variant 
        overallScore
        distances {
          sourceId
          aggregatedScore
          tissues {
            distance
          }
        }
      }
    }"""
    # Set base URL of GraphQL API endpoint
    base_url = "https://api.genetics.opentargets.org/graphql"
    # Set variables object of arguments to be passed to endpoint
    variables = {'variantId': variantId}
    # Construct POST request body object with query string and variables
    json={"query": query_string, "variables": variables}

    # Perform POST request
    r = requests.post(url=base_url, json=json)
    print(r.status_code)
    r.json()

    return r 

st.cache_data
def transform_json_to_df_open_targets(r):
    """Turning output of OpenTargets VSG into a dataframe. """

    tmp = pd.DataFrame(r.json()["data"]["genesForVariant"])
    tmp["gene"] = [i["id"] for i in tmp["gene"]]

    tmp["distances"] = [i[0] if len(i)>0 else {} for i in tmp["distances"]]
    tissues_df = pd.json_normalize(tmp["distances"])
    tissues_df["tissues"] = tissues_df["tissues"].apply(lambda x: {'distance': x[0]['distance']} if isinstance(x, list) and len(x) > 0 else {})

    distances_df = pd.json_normalize(tissues_df["tissues"])
    tissues_df = pd.merge(tissues_df, distances_df, right_index=True, left_index=True )
    tissues_df = tissues_df.drop(['sourceId', 'tissues'], axis=1)
    
    qtls = [i if i else [] for i in tmp["qtls"]]

    normalized_data = []
    for item in qtls:
        if not item:
            normalized_data.append({})
        else:
            normalized_data.append({entry['typeId']: entry['aggregatedScore'] for entry in item})

    # Create a DataFrame
    normalised_qtls = pd.DataFrame(normalized_data)

    tmp =pd.merge(tmp, tissues_df, left_index=True, right_index=True)
    tmp = pd.merge(tmp, normalised_qtls, left_index=True, right_index=True)
    tmp = tmp.drop(["qtls", "distances"], axis =1)
    return tmp

def open_targets_df(variant_code_b38_ref_alt):
    r = open_targets_v2g(variant_code_b38_ref_alt)
    df = transform_json_to_df_open_targets(r)
    return df


st.cache_data
# Open Targets PheWAS
def OT_phewas(variant_code_b38_ref_alt):
    query_string = """
    query PheWASQuery($variantId: String!) {
        pheWAS(variantId: $variantId){
        associations {
            study {
                traitReported
                traitCategory
                pmid
                pubDate
                pubAuthor
                source
                }
            pval
            beta
            oddsRatio
            nTotal
            nCases
            eaf
            se
            }
        }
    }"""

    # Set base URL of GraphQL API endpoint
    base_url = "https://api.genetics.opentargets.org/graphql"
    variantId = variant_code_b38_ref_alt
    # Set variables object of arguments to be passed to endpoint
    variables = {'variantId': variantId}
    # Construct POST request body object with query string and variables
    json={"query": query_string, "variables": variables}

    # Perform POST request
    r = requests.post(url=base_url, json=json)
    print(r.status_code)
    r.json()
    return r
st.cache_data
def unravel_phewas_open_targets(r):
    tmp = pd.DataFrame(r.json()["data"]["pheWAS"])
    tmp["study"] = [i["study"] for i in tmp["associations"]]
    tmp["pval"] =[i["pval"] for i in tmp["associations"]]
    tmp["oddsRatio"] =[i["oddsRatio"] for i in tmp["associations"]]
    tmp["nTotal"] =[i["nTotal"] for i in tmp["associations"]]
    tmp["nCases"] =[i["nCases"] for i in tmp["associations"]]
    tmp["eaf"] =[i["eaf"] for i in tmp["associations"]]
    tmp["se"] =[i["se"] for i in tmp["associations"]]
    tmp["traitReported"] = [i["traitReported"] for i in tmp["study"]]
    tmp["pmid"] = [i["pmid"] for i in tmp["study"]]
    tmp["pubDate"] = [i["pubDate"] for i in tmp["study"]]
    tmp["pubAuthor"] = [i["pubAuthor"] for i in tmp["study"]]
    tmp["source"]= [i["source"] for i in tmp["study"]]
    tmp["traitCategory"] = [i["traitCategory"] for i in tmp["study"]]
    tmp = tmp.drop(columns=['study', 'associations'])
    return tmp

def phewas(variant_code_b38_ref_alt):
    r = OT_phewas(variant_code_b38_ref_alt)
    df = unravel_phewas_open_targets(r)
    return df

