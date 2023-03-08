import streamlit as st
import plotly.express as px
import pandas as pd
import numpy as np
from tarviz.utils import load_data, bQTLs_data
from tarviz.widgets import pvalue_filters_widget, filter, top_page_widget
from argparse import ArgumentParser

st.set_page_config(layout="wide", page_icon="images/logo.ico")

@st.cache_data
def target_options(data):
    return ["None", *data.TARGET.unique()]

def treatment_options(data):
    return ["None", *data.TREATMENTS.unique()]

@st.cache_data
def pvalues_hist(data: pd.DataFrame):
    data["NEG_LOG_PVALUE"] = - np.log(data["PVALUE"])
    fig = px.histogram(
        data, 
        x="NEG_LOG_PVALUE", 
        labels={"NEG_LOG_PVALUE": f"-log(PVALUE)"},
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

    data = load_data(bQTLs_data())

    # Title
    top_page_widget()

    st.markdown("Welcome to the TarGene visualization interface.")
    # Display the table
    mt_method, pvalue = pvalue_filters_widget()
    col21, colr22 = st.columns([2, 1])
    col211, col212 = col21.columns(2)
    target = col211.selectbox("Target", target_options(data))
    treatment_combo = col212.selectbox("Treatment", treatment_options(data))
    # Data Table
    filtered = filter(data, mt_method, pvalue, target, treatment_combo, "None")
    col21.markdown("Only first 1000 results are presented.")
    col21.dataframe(limit_data(filtered).style.hide(axis='index'), use_container_width=True)
    # Hits per Target/Treatment
    per = colr22.multiselect("Hit counts", ["TREATMENTS", "TARGET"], default=["TREATMENTS"])
    nhits = filtered.groupby(per).size().reset_index(name='COUNTS').sort_values("COUNTS", ascending=False)
    colr22.dataframe(nhits.style.hide(axis='index'))
    # Plot Pvalues
    pvalues_hist(filtered)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("rundir", type=str, default=".", help='Directory where the targene-pipeline was run from.')
    args = parser.parse_args()
    main(args)

