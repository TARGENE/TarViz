import streamlit as st
import plotly.express as px
import pandas as pd
import numpy as np
from tarviz.utils import load_data, bQTLs_data
from tarviz.widgets import sidebar_widget, filter, top_page_widget
from argparse import ArgumentParser

st.set_page_config(layout="wide", page_icon="images/logo.ico")

@st.cache_data
def target_options(data):
    return ["None", *data.TARGET.unique()]

def treatment_options(data):
    return ["None", *data.TREATMENTS.unique()]

@st.cache_data
def pvalues_hist(data: pd.DataFrame, scale, pval_col):
    data[f"NEG_LOG_{pval_col}"] = - np.log(data[pval_col])
    if scale == "log":
        x = f"NEG_LOG_{pval_col}"
        label = "-log(f{pval_col})"
    else:
        x = pval_col
        label = pval_col
    fig = px.histogram(
        data, 
        x=x, 
        labels={x: label},
        nbins=100
    )
    st.plotly_chart(fig, use_container_width=True)

@st.cache_data
def limit_data(data, maxrows=1000):
    return data[:maxrows]

@st.cache_data
def columns(data):
    return data.columns

def main(args):
    st.session_state['nextflow_rundir'] = args.rundir

    results_file, pval_col, pvalue = sidebar_widget()

    data = load_data(bQTLs_data(), results_file)

    # Title
    top_page_widget()

    st.markdown("Welcome to the TarGene visualization interface.")
    # Display the table
    col21, colr22 = st.columns([2, 1])
    col211, col212 = col21.columns(2)

    target_from_selectbox = col211.selectbox("Select a target", target_options(data))
    target_from_text = col211.text_input("or provide a regular expression", "None", key=1)
    target_filter = target_from_selectbox if target_from_selectbox != "None" else target_from_text

    treatment_from_selectbox = col212.selectbox("Treatment", treatment_options(data))
    treatment_from_text = col212.text_input("or provide a regular expression", "None", key=2)
    treatment_filter = treatment_from_selectbox if treatment_from_selectbox != "None" else treatment_from_text

    # Data Table
    filtered = filter(data, pval_col, pvalue, target_filter, treatment_filter, "None")
    col21.markdown("Only first 1000 results are presented.")
    col21.dataframe(limit_data(filtered).style.hide(axis='index'), use_container_width=True)
    # Hits per Target/Treatment
    per = colr22.multiselect("Hit counts", ["TREATMENTS", "TARGET"], default=["TREATMENTS"])
    nhits = filtered.groupby(per).size().reset_index(name='COUNTS').sort_values("COUNTS", ascending=False)
    colr22.dataframe(nhits.style.hide(axis='index'))
    # Plot Pvalues
    scale = st.selectbox("Scale", ["None", "log"])
    pvalues_hist(filtered, scale, pval_col)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("rundir", type=str, default=".", help='Directory where the targene-pipeline was run from.')
    args = parser.parse_args()
    main(args)