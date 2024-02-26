import streamlit as st
from st_aggrid import AgGrid, GridOptionsBuilder
from tarviz.utils import bQTLs_data, http_ensembl_binding_matrix, bqtls_hit_counts, \
    load_data, feature_columns, look_up_variant_gtex, open_targets_df, phewas, look_up_variant_gtex_tissue
from tarviz.widgets import sidebar_widget, filter, \
    SNPinfo, region_annotations, top_page_widget, plot_motif_logo, SNPinfo_priority
from tarviz.constants import ANNOTATION_FEATURES, PRIORITY_FEATURES
import tables
import pandas as pd

#modulation_plot,
top_page_widget()
results_file, mt_method, pvalue = sidebar_widget()
bqtls_data = bQTLs_data()

bqtl_df = pd.DataFrame(bqtls_data)
st.dataframe(bqtl_df)
bqtl = bqtl_df["ID"].iloc[1]
bqtls = bqtl_df["ID"].tolist()[0:19]


data = filter(load_data(bqtls_data, results_file), mt_method, pvalue, "None", "None", "None")


def no_index(df):
    return df.reset_index(drop=True)

st.subheader("Priority annotations")

col1, col2 = st.columns(2)
distance = col1.number_input("Region size (bp)", min_value=0, max_value=5000000, value=10)
v2g = col1.number_input("V2G score ", min_value=0.0, max_value=1.0, value=0.5)
d2g = col1.number_input("Distance to gene ", min_value=0, max_value=5000000, value=5000)
features = col2.multiselect("Features", ANNOTATION_FEATURES, default=["gene", "regulatory"])
priority = col2.multiselect("Extra priorities", PRIORITY_FEATURES, default=["trait", "V2G"])
trait = col2.text_input("Trait of interest", "None")

annotated_list =[]
non_annotated_list =[]

for b in bqtls:
    for feature in features:
        basesnpinfo = SNPinfo_priority(b, bqtls_data)
        chr = basesnpinfo["Chromosome"][0]
        v_start = basesnpinfo["Start"][0]
        v_end = basesnpinfo["End"][0]
        variant_code_b38_ref_alt = str(basesnpinfo["Chromosome"][0])+"_"+str(basesnpinfo["Start"][0])+"_"+str(basesnpinfo["REF Allele"][0])+"_"+str(basesnpinfo["ALT Allele"][0])
        annotations = region_annotations(chr, v_start, v_end, distance, (feature,))



        if len(annotations) > 0:
            annotations = feature_columns(annotations, feature)
            
            open_targets = open_targets_df(variant_code_b38_ref_alt)
            ot_phewas= phewas(variant_code_b38_ref_alt)
           


            #### Selection
            selected_open_targets = open_targets[(open_targets["aggregatedScore"] > v2g)&(open_targets["distance"] < d2g)]


            selected_ot_phewas = ot_phewas[ot_phewas["traitReported"].str.contains(trait, case=False)]
            basesnpinfo_replicated = pd.concat([basesnpinfo] * len(annotations), ignore_index=True)
            ot_phews_base = pd.concat([no_index(basesnpinfo_replicated),no_index(selected_ot_phewas)],axis =1)



            df = pd.concat([no_index(basesnpinfo_replicated),
                            no_index(annotations),
                           no_index(selected_open_targets)], 
                         #  no_index(selected_ot_phewas)],
                          axis=1)
            final_df_annotated = df.merge(ot_phews_base , how='inner', on=basesnpinfo.columns.tolist())
    
            annotated_list.append(final_df_annotated)
 
    
        


  
        else:
            open_targets = open_targets_df(variant_code_b38_ref_alt)
            ot_phewas= phewas(variant_code_b38_ref_alt)
           


            #### Selection
            selected_open_targets = open_targets[(open_targets["aggregatedScore"] > v2g)&(open_targets["distance"] < d2g)]


            selected_ot_phewas = ot_phewas[ot_phewas["traitReported"].str.contains(trait, case=False)]
            #basesnpinfo_replicated = pd.concat([basesnpinfo] * len(annotations), ignore_index=True)
            ot_phews_base = pd.concat([no_index(basesnpinfo),no_index(selected_ot_phewas)],axis =1)



            df = pd.concat([no_index(basesnpinfo),
                           no_index(selected_open_targets)], 
                         #  no_index(selected_ot_phewas)],
                          axis=1)
            final_df_non = df.merge(ot_phews_base , how='inner', on=basesnpinfo.columns.tolist())
            non_annotated_list.append(final_df_non)
    


non_anno_df = pd.concat(non_annotated_list, axis =0)
anno_df = pd.concat(annotated_list, axis =0)

priority_df = pd.concat([anno_df, non_anno_df] ,axis =0)
st.dataframe(priority_df)

            
    
       




