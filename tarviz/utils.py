import pandas as pd
import os
import requests
import streamlit as st
from tarviz.constants import ENSEMBL_URL

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

@st.cache_data
def bQTLs_data(nextflow_rundir, config_filename="nextflow.config"):
    params = load_pipeline_params(os.path.join(nextflow_rundir, config_filename))
    bqtls_filepath = os.path.join(nextflow_rundir, params["BQTLS"])
    return pd.read_csv(bqtls_filepath)

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

