import streamlit as st
from st_aggrid import AgGrid, GridOptionsBuilder
from tarviz.utils import bQTL_list, raw_data_file
from tarviz.widgets import pvalue_filters_widget, filter, load_data, modulation_plot, SNPinfo

st.set_page_config(layout="wide", page_icon="images/logo.ico")

st.title("bQTL View")
bqtl = st.selectbox("bQTL rsID: ", bQTL_list(st.session_state.nextflow_rundir))
mt_method, pvalue = pvalue_filters_widget()

SNPinfo(bqtl)

st.header(f"TarGene Hits")
st.markdown("Please select an estimation result to display further information (this may take a few seconds):")
data = filter(load_data(), mt_method, pvalue, "None", "None", bqtl)

builder = GridOptionsBuilder.from_dataframe(data)
builder.configure_selection(selection_mode='single', use_checkbox=True, pre_selected_rows=[0])
builder.configure_pagination(enabled=True, paginationAutoPageSize=True, paginationPageSize=10)
grid_return = AgGrid(data, gridOptions=builder.build())

if len(grid_return["selected_rows"]) > 0:
    modulation_plot(bqtl, grid_return["selected_rows"][0])