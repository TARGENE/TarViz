import streamlit as st
from tarviz.utils import bQTL_list

st.set_page_config(layout="wide", page_icon="images/logo.ico")
NEXTFLOW_RUNDIR = "tests/nextflow_rundir"

st.selectbox("bQTL rsID: ", bQTL_list(NEXTFLOW_RUNDIR))