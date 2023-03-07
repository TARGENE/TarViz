import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from tarviz.constants import MT_ADJUSTEMENT_METHODS, DATA_COLUMNS
from tarviz.utils import result_file, raw_data_file, http_variant_info, http_ensemble_annotations


def top_page_widget():
    col1, col2 = st.columns(2)
    col1.title("TarGene Dashboard")
    col2.image("images/logo.jpg")

def pvalue_filters_widget():
    mt_method = st.sidebar.selectbox("Multiple Testing Adjustment Method", MT_ADJUSTEMENT_METHODS)
    pvalue = float(st.sidebar.text_input("Pvalue Threshold", value=0.05))
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

@st.cache_data
def modulation_plot(bqtl, selected):
    trait = selected["TARGET"]
    # Is IATE
    if selected["PARAMETER_TYPE"] == "IATE":
        variants = selected["TREATMENTS"].split("_&_")
        eqtl = variants[1] if variants[0] == bqtl else variants[0]
        v1_case, v2_case = (float(x) for x in selected["CASE"].split("_&_"))
        v1_control, v2_control = (float(x) for x in selected["CONTROL"].split("_&_"))

    raw_data = pd.read_csv(
        raw_data_file(st.session_state.nextflow_rundir), 
        usecols = [trait, *(x for x in variants)],
    ).dropna()
    raw_data[bqtl] = raw_data[bqtl].astype(int)
    raw_data[eqtl] = raw_data[eqtl].astype(int)
    raw_data = raw_data[(raw_data[variants[0]].isin([v1_case, v1_control])) & (raw_data[variants[1]].isin([v2_case, v2_control]))]

    if raw_data[trait].nunique() == 2:
        means = raw_data.groupby(variants).agg(
            mean=(trait, np.mean),
            ncases=(trait, np.sum)
            ).reset_index()
        means[eqtl] = means[eqtl].astype(str)
        fig = px.scatter(
            means, 
            y="mean", 
            x=bqtl, 
            color=eqtl, 
            size="ncases", 
            labels={
                    "mean": f"{trait} mean",
                    bqtl: f"{bqtl} (# minor alleles)",
                    eqtl: f"{eqtl} (# minor alleles)"
            },
            title=f"{bqtl}'s modulation by {eqtl} (from raw data)",
            color_discrete_sequence=["red", "blue"]
        )
        fig.update_layout(
            scattermode="group", 
            scattergap=0.9,
            xaxis_type='category',
            title_x=0.35
        )
        
    # Continuous traits
    else:
        fig = px.violin(
            raw_data, 
            y=trait, 
            x=bqtl, 
            color=eqtl,
            labels={
                    bqtl: f"{bqtl} (# minor alleles)",
                    eqtl: f"{eqtl} (# minor alleles)"
            }, 
            title=f"{bqtl}'s modulation by {eqtl} (from raw data)",
            box=True,
            color_discrete_sequence=["red", "blue"]
            )
        fig.update_layout(
            violinmode="group", 
            violingap=0.7,
            xaxis_type='category',
            title_x=0.35
        )

    st.plotly_chart(fig, use_container_width=True)

@st.cache_data
def SNPinfo(rsid):
    response = http_variant_info(rsid)
    mapping_1 = response["mappings"][0]
    basesnpinfo = {
        "Location": [mapping_1["location"]],
        "Strand": [mapping_1["strand"]],
        "Alleles": [mapping_1["allele_string"]],
        "Minor Allele": [response["minor_allele"]],
        "MAF": [response["MAF"]],
        "Ancestral Allele": [mapping_1["ancestral_allele"]],
    }
    st.header(f"Ensembl Report")
    # Inject CSS with Markdown to hide table's row indices
    hide_table_row_index = """
            <style>
            thead tr th:first-child {display:none}
            tbody th {display:none}
            </style>
            """
    st.markdown(hide_table_row_index, unsafe_allow_html=True)
    st.table(basesnpinfo)
    st.subheader("Existing phenotyping annotations")
    if len(response["phenotypes"]) > 0:
        st.dataframe(response["phenotypes"])
    else:
        st.markdown("None.")
    
    return basesnpinfo

@st.cache_data
def region_annotations(location_str, distance, features):
    response = http_ensemble_annotations(
        location_str, 
        distance=distance, 
        features=features
        )
    return pd.DataFrame(response)
    
