import pandas as pd
import numpy as np
import streamlit as st
import plotly.express as px
from st_aggrid import AgGrid, GridOptionsBuilder
from tarviz.utils import bQTL_list, raw_data_file
from tarviz.widgets import pvalue_filters_widget, filter, load_data

st.set_page_config(layout="wide", page_icon="images/logo.ico")

st.title("bQTL View")
bqtl = st.selectbox("bQTL rsID: ", bQTL_list(st.session_state.nextflow_rundir))
mt_method, pvalue = pvalue_filters_widget()

data = filter(load_data(), mt_method, pvalue, "None", "None", bqtl)

builder = GridOptionsBuilder.from_dataframe(data)
builder.configure_selection(selection_mode='single', use_checkbox=True, pre_selected_rows=[0])
builder.configure_pagination(enabled=True, paginationAutoPageSize=True, paginationPageSize=20)
grid_return = AgGrid(data, gridOptions=builder.build())

selected = grid_return["selected_rows"][0]
trait = selected["TARGET"]
# Is IATE
if selected["PARAMETER_TYPE"] == "IATE":
    variants = selected["TREATMENTS"].split("_&_")
    eqtl = variants[1] if variants[0] == bqtl else variants[0]
    cols_to_load = [trait, *(v for v in variants)]
    v1_case, v2_case = (float(x) for x in selected["CASE"].split("_&_"))
    v1_control, v2_control = (float(x) for x in selected["CONTROL"].split("_&_"))

raw_data = pd.read_csv(
        raw_data_file(st.session_state.nextflow_rundir), 
        usecols = [trait, *(x for x in variants)],
    ).dropna()
raw_data[bqtl] = raw_data[bqtl].astype(int)
raw_data[eqtl] = raw_data[eqtl].astype(int)
raw_data = raw_data[(raw_data[variants[0]].isin([v1_case, v1_control])) & (raw_data[variants[1]].isin([v2_case, v2_control]))]

# Binary traits
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