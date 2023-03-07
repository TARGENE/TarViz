import streamlit as st
from st_aggrid import AgGrid, GridOptionsBuilder
from tarviz.utils import bQTL_list, raw_data_file
from tarviz.widgets import pvalue_filters_widget, filter, load_data, modulation_plot, \
    SNPinfo, region_annotations, top_page_widget
from tarviz.constants import ANNOTATION_FEATURES


top_page_widget()

bqtl = st.selectbox("bQTL", bQTL_list(st.session_state.nextflow_rundir))
mt_method, pvalue = pvalue_filters_widget()

# Display SNP base info
basesnpinfo = SNPinfo(bqtl)
# Display region information
st.subheader("Region annotations")
col1, col2 = st.columns(2)
distance = col1.number_input("Region size (bp)", min_value=0, max_value=5000000, value=10)
features = col2.multiselect("Features", ANNOTATION_FEATURES, default=["gene", "regulatory", "motif"])
tabs = st.tabs(features)
for (tab, feature) in zip(tabs, features):
    location = basesnpinfo["Location"][0]
    annotations = region_annotations(location, distance, (feature,))
    with tab:
        if len(annotations) > 0:
            if feature == "motif":
                _, v_loc = location.split(":")
                v_min_loc, v_max_loc = (int(x) for x in v_loc.split("-"))
                highlighted_rows = ((annotations['start'] <= v_min_loc) & (annotations['end'] >= v_max_loc)).map({
                    True: 'background-color: green',
                    False: ''
                })
                st.markdown("Highlighted rows contain the variant.")
                st.dataframe(annotations.style.apply(lambda _: highlighted_rows))
                st.selectbox("Binding Matrix", annotations["binding_matrix_stable_id"].unique())
            else:
                st.dataframe(annotations)
        else:
            st.markdown("No annotation found.")

st.header(f"TarGene Hits")
st.markdown("Please select an estimation result to display further information (this may take a few seconds):")
data = filter(load_data(), mt_method, pvalue, "None", "None", bqtl)

builder = GridOptionsBuilder.from_dataframe(data)
builder.configure_selection(selection_mode='single', use_checkbox=True, pre_selected_rows=[0])
builder.configure_pagination(enabled=True, paginationAutoPageSize=False, paginationPageSize=10)
grid_return = AgGrid(data, gridOptions=builder.build())

if len(grid_return["selected_rows"]) > 0:
    modulation_plot(bqtl, grid_return["selected_rows"][0])