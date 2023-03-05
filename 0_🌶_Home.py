import streamlit as st
import plotly.express as px
import pandas as pd
import base64

st.set_page_config(layout="wide", page_icon="images/logo.ico")

COLUMNS = [
    'PARAMETER_TYPE', 'TREATMENTS', 'TARGET', 'PVALUE', 'ADJUSTED_PVALUE', 
    'ESTIMATE', 'LWB', 'UPB', 'STD', 'CASE', 'CONTROL', 'CONFOUNDERS', 
    'COVARIATES', 'INITIAL_ESTIMATE', 'LOG']

MT_ADJUSTEMENT_METHODS = ["Trait based Benjamini Hochberg", "None"]

@st.cache_data
def target_options(data):
    return ["None", *data.TARGET.unique()]

def treatment_options(data):
    return ["None", *data.TREATMENTS.unique()]

@st.cache_data
def load_data():
    return pd.read_csv("preliminary_results.csv")[COLUMNS]

@st.cache_data
def filter(df: pd.DataFrame, mt_method, pvalue, target, treatment) -> pd.DataFrame:
    # Pvalue based filter
    filterstring = f"(df.ADJUSTED_PVALUE < {pvalue})"
    if mt_method == "None":
        filterstring = f"(df.PVALUE < {pvalue})"
    # Target based filter
    if target != "None":
        filterstring += f" & (df.TARGET == \"{target}\")"
    # Treatment based filter
    if treatment != "None":
        filterstring += f" & (df.TREATMENTS == \"{treatment}\")"
    
    return pd.eval(f"df[{filterstring}]")

@st.cache_data
def pvalues_hist(data: pd.DataFrame):
    fig = px.histogram(data, x="PVALUE", nbins=100)
    st.plotly_chart(fig, use_container_width=True)

@st.cache_data
def limit_data(data, maxrows=1000):
    return data[:maxrows]

@st.cache_data
def columns(data):
    return data.columns

data = load_data()

# Title
st.title("TarGene Dashboard")
st.markdown("Welcome to the TarGene visualization interface.")
# Display the table
col1, col2 = st.columns(2)
mt_method = col1.selectbox("Multiple Testing Adjustment Method", MT_ADJUSTEMENT_METHODS)
pvalue = float(col2.text_input("Pvalue Threshold", value=0.05))
st.markdown("""---""")
col21, colr22 = st.columns([2, 1])
col211, col212 = col21.columns(2)
target = col211.selectbox("Target", target_options(data))
treatment = col212.selectbox("Treatment", treatment_options(data))
# Data Table
filtered = filter(data, mt_method, pvalue, target, treatment)
col21.markdown("Only first 1000 results are presented.")
col21.dataframe(limit_data(filtered).style.hide(axis='index'), use_container_width=True)
# Hits per Target/Treatment
per = colr22.multiselect("Hit counts", ["TREATMENTS", "TARGET"], default=["TREATMENTS"])
print(per)
nhits = filtered.groupby(per).size().reset_index(name='COUNTS')
colr22.dataframe(nhits.style.hide(axis='index'))
# Plot Pvalues
pvalues_hist(filtered)

