import pandas as pd
import os
import requests
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
    data = pd.read_csv(os.path.join(st.session_state["nextflow_rundir"], "results", results_file))[DATA_COLUMNS]
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
    st.write(counts)
    return counts["BQTL"] + " (" + counts["COUNTS"].astype(str) + " hits)"

@st.cache_data
def feature_columns(df, feature):
    if feature == "motif":
        return df[["transcription_factor_complex", "binding_matrix_stable_id", "score", "start", "end", "strand"]]
    else:
        return df





