import pandas as pd
import os

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
    
def bQTL_list(nextflow_rundir, config_filename="nextflow.config"):
    params = load_pipeline_params(os.path.join(nextflow_rundir, config_filename))
    bqtls_filepath = os.path.join(nextflow_rundir, params["BQTLS"])
    return pd.read_csv(bqtls_filepath, usecols=["ID"]).ID.tolist()
