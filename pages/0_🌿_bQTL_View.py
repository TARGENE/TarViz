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
# Is ATE
else:
    variants = [bqtl]

raw_data = pd.read_csv(
        raw_data_file(st.session_state.nextflow_rundir), 
        usecols = [trait, *(x for x in variants)]
    )

gpd = raw_data.groupby(variants)
# Binary traits
if raw_data[trait].nunique() == 2:
    means = gpd.agg(
        Mean=(trait, np.mean),
        Sum=(trait, np.sum)
        ).reset_index()
    means[eqtl] = means[eqtl].astype(str)
    means["LogSum"] = np.log(means["Sum"] + 1)
    st.dataframe(means)
    fig = px.scatter(means, y="Mean", x=bqtl, color=eqtl, size="Sum")
    fig.update_layout(scattermode="group", scattergap=0.9)
    
# Continuous traits
else:
    gpd = raw_data.groupby([bqtl])

st.plotly_chart(fig, use_container_width=True)