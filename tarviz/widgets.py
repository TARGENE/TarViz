import streamlit as st
import pandas as pd
from tarviz.constants import MT_ADJUSTEMENT_METHODS, DATA_COLUMNS

from tarviz.utils import result_file

def pvalue_filters_widget():
    col1, col2 = st.columns(2)
    mt_method = col1.selectbox("Multiple Testing Adjustment Method", MT_ADJUSTEMENT_METHODS)
    pvalue = float(col2.text_input("Pvalue Threshold", value=0.05))
    st.markdown("""---""")
    return mt_method, pvalue

@st.cache_data
def load_data():
    return pd.read_csv(result_file(st.session_state.nextflow_rundir))[DATA_COLUMNS]

@st.cache_data
def unique_treatments(data):
    st.selectbox("Treatment", [t for x in data.TREATMENTS for t in x.split("_&_")])

@st.cache_data
def filter(df: pd.DataFrame, mt_method, pvalue, target, treatment_combo, treatment) -> pd.DataFrame:
    # Pvalue based filter
    filterstring = f"(df.ADJUSTED_PVALUE < {pvalue})"
    if mt_method == "None":
        filterstring = f"(df.PVALUE < {pvalue})"
    # Target based filter
    if target != "None":
        filterstring += f" & (df.TARGET == \"{target}\")"
    # Treatment combo based filter
    if treatment_combo != "None":
        filterstring += f" & (df.TREATMENTS == \"{treatment_combo}\")"
    # Treatment based filter
    if treatment != "None":
        filterstring += f" & (df.TREATMENTS.str.contains(\"{treatment}\"))"
    
    return pd.eval(f"df[{filterstring}]")