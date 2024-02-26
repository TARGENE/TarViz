import streamlit as st
from st_aggrid import AgGrid, GridOptionsBuilder
from tarviz.utils import bQTLs_data, http_ensembl_binding_matrix, bqtls_hit_counts, \
    load_data, feature_columns, look_up_variant_gtex, open_targets_df, phewas, look_up_variant_gtex_tissue
from tarviz.widgets import sidebar_widget, filter, \
    SNPinfo, region_annotations, top_page_widget, plot_motif_logo
from tarviz.constants import ANNOTATION_FEATURES
import tables

#modulation_plot,
top_page_widget()
results_file, mt_method, pvalue = sidebar_widget()
bqtls_data = bQTLs_data()

data = filter(load_data(bqtls_data, results_file), mt_method, pvalue, "None", "None", "None")


bqtl_str = st.selectbox("Select bQTL (# Trait hits)", bqtls_hit_counts(data))
bqtl = bqtl_str.split(" (")[0]
# Display SNP base info
basesnpinfo = SNPinfo(bqtl, bqtls_data)
# Display region information
st.subheader("Region annotations")
col1, col2 = st.columns(2)
distance = col1.number_input("Region size (bp)", min_value=0, max_value=5000000, value=10)
features = col2.multiselect("Features", ANNOTATION_FEATURES, default=["gene", "regulatory", "motif"])
tabs = st.tabs(features)
for (tab, feature) in zip(tabs, features):
    chr = basesnpinfo["Chromosome"][0]
    v_start = basesnpinfo["Start"][0]
    v_end = basesnpinfo["End"][0]
    variant_code_b38_ref_alt = str(basesnpinfo["Chromosome"][0])+"_"+str(basesnpinfo["Start"][0])+"_"+str(basesnpinfo["REF Allele"][0])+"_"+str(basesnpinfo["ALT Allele"][0])
    annotations = region_annotations(chr, v_start, v_end, distance, (feature,))

    with tab:
        if len(annotations) > 0:
            annotations = feature_columns(annotations, feature)
            st.markdown("Highlighted rows correspond to regions containing the variant.")
            highlighted_rows = ((annotations['start'] <= v_start) & (annotations['end'] >= v_end)).map({
                        True: 'background-color: green',
                        False: ''
            })
            st.dataframe(annotations.style.apply(lambda _: highlighted_rows))
            if feature == "motif":
                matrix_id = st.selectbox("Binding Matrix", annotations["binding_matrix_stable_id"].unique())
                response = http_ensembl_binding_matrix(matrix_id)                
                binding_motif_start = annotations[annotations.binding_matrix_stable_id == matrix_id].iloc[0]["start"]
                pmin = int(v_start - binding_motif_start + 1)
                pmax = int(v_end - binding_motif_start + 1)
                plot_motif_logo(response, pmin, pmax)
            



        else:
            st.markdown("No annotation found.")

# OpenTargets V2G

st.header("Open Targets")
st.subheader("Open Targets V2G")
st.dataframe(open_targets_df(variant_code_b38_ref_alt))


st.subheader("Open Targets PheWas")
st.dataframe(phewas(variant_code_b38_ref_alt))

st.header("GTEx expression")
gtex, selection_tissue = look_up_variant_gtex_tissue(variant_code_b38_ref_alt) 

selection_tissue = st.selectbox("Select tissue",selection_tissue)
st.dataframe(look_up_variant_gtex(gtex, variant_code_b38_ref_alt, selection_tissue))


st.header(f"TarGene Hits")
st.markdown("Please select an estimation result to display further information (this may take a few seconds):")

bqtl_data = filter(data, mt_method, pvalue, "None", "None", bqtl)
builder = GridOptionsBuilder.from_dataframe(bqtl_data)
builder.configure_selection(selection_mode='single', use_checkbox=True, pre_selected_rows=[0])
builder.configure_pagination(enabled=True, paginationAutoPageSize=False, paginationPageSize=10)
grid_return = AgGrid(bqtl_data, gridOptions=builder.build())

#if len(grid_return["selected_rows"]) > 0:
#    modulation_plot(bqtl, grid_return["selected_rows"][0])