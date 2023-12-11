import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import logomaker
from tarviz.constants import MT_ADJUSTEMENT_METHODS
from tarviz.utils import raw_data_file, http_variant_info, http_ensemble_annotations, look_up_variant_gtex

def normal_approx_error(m, n):
    return 1.96 * np.sqrt(np.divide(m * (1 - m), n))

def top_page_widget():
    col1, col2 = st.columns(2)
    col1.title("TarGene Dashboard")
    col2.image("images/logo.jpg")

def sidebar_widget():
    results_file = st.sidebar.text_input("Result's filename", value="summary.csv")
    column = st.sidebar.selectbox("P-value column", MT_ADJUSTEMENT_METHODS)
    pvalue = float(st.sidebar.text_input("Pvalue Threshold", value=0.05))
    return results_file, column, pvalue

@st.cache_data
def unique_treatments(data):
    st.selectbox("Treatment", [t for x in data.TREATMENTS for t in x.split("_&_")])

@st.cache_data
def filter(df: pd.DataFrame, pval_col, pvalue, target, treatment_combo, treatment) -> pd.DataFrame:
    # Pvalue based filter
    filterstring = f"({pval_col} < {pvalue})"
    # Target based filter
    if target != "None":
        filterstring += f" & (TARGET.str.contains(\"{target}\"))"
    # Treatment combo based filter
    if treatment_combo != "None":
        filterstring += f" & (TREATMENTS.str.contains(\"{treatment_combo}\"))"
    # Treatment based filter
    if treatment != "None":
        filterstring += f" & (TREATMENTS.str.contains(\"{treatment}\"))"

    st.write(filterstring)

    return df.query(filterstring)

def binary_IATE_plot(raw_data, bqtl, modulator, trait):
    means = raw_data.groupby([bqtl, modulator]).agg(
        mean=(trait, np.mean),
        count=(trait, np.size),
        ncases=(trait, np.sum)
    ).reset_index()
    means["Error"] = normal_approx_error(means["mean"], means["count"])
    means[modulator] = means[modulator].astype(str)
    fig = px.scatter(
        means, 
        error_y="Error",
        y="mean", 
        x=bqtl, 
        color=modulator, 
        size="ncases", 
        labels={
                "mean": f"{trait} mean",
                bqtl: f"{bqtl} (# minor alleles)",
                modulator: f"{modulator} (# minor alleles)"
        },
        title=f"{bqtl}'s modulation by {modulator} (from raw data).<br><sup>Error bars are based on the Normal approximation.</sup>",
        color_discrete_sequence=["red", "blue"]
        )
    fig.update_layout(
            scattermode="group", 
            scattergap=0.9,
            xaxis_type='category',
            title_x=0.35
        )
    return fig

def continuous_IATE_plot(raw_data, bqtl, modulator, trait):
    fig = px.violin(
        raw_data, 
        y=trait, 
        x=bqtl, 
        color=modulator,
        labels={
                bqtl: f"{bqtl} (# minor alleles)",
                modulator: f"{modulator} (# minor alleles)"
        }, 
        title=f"{bqtl}'s modulation by {modulator} (from raw data)",
        box=True,
        color_discrete_sequence=["red", "blue"]
    )
    fig.update_layout(
        violinmode="group", 
        violingap=0.7,
        xaxis_type='category',
        title_x=0.35
    )
    return fig

def binary_ATE_plot(raw_data, bqtl, trait):
    means = raw_data.groupby([bqtl]).agg(
        mean=(trait, np.mean),
        count=(trait, np.size),
        ncases=(trait, np.sum)
    ).reset_index()
    means["Error"] = normal_approx_error(means["mean"], means["count"])
    fig = px.scatter(
        means, 
        y="mean",
        error_y="Error",
        x=bqtl, 
        size="ncases", 
        labels={
                "mean": f"{trait} mean",
                bqtl: f"{bqtl} (# minor alleles)",
        },
        title=f"{bqtl}'s effect (from raw data).<br><sup>Error bars are based on the Normal approximation.</sup>",
    )
    fig.update_layout(
        scattermode="group", 
        scattergap=0.9,
        xaxis_type='category',
        title_x=0.35
    )
    return fig

def continuous_ATE_plot(raw_data, bqtl, trait):
    fig = px.violin(
        raw_data, 
        y=trait, 
        x=bqtl, 
        labels={
                bqtl: f"{bqtl} (# minor alleles)",
        }, 
        title=f"{bqtl}'s effect (from raw data)",
        box=True,
    )
    fig.update_layout(
        violinmode="group", 
        violingap=0.7,
        xaxis_type='category',
        title_x=0.35
    )
    return fig

"""@st.cache_data
def modulation_plot(bqtl, selected):
    trait = selected["TARGET"]
    param_type = selected["PARAMETER_TYPE"]
    treatments_string = selected["TREATMENTS"]

    if param_type == "IATE":
        variants = treatments_string.split("_&_")
        if len(variants) != 2:
            st.markdown("The modulation plot for the IATE is currently only available when the Treatment contains exactly 2 variables.")
            return
        else:
            modulator = variants[1] if variants[0] == bqtl else variants[0]
            v1_case, v2_case = (float(x) for x in selected["CASE"].split("_&_"))
            v1_control, v2_control = (float(x) for x in selected["CONTROL"].split("_&_"))
            raw_data = pd.read_csv(
                raw_data_file(st.session_state.nextflow_rundir), 
                usecols = [trait, *(x for x in variants)],
            ).dropna()
            raw_data[bqtl] = raw_data[bqtl].astype(int)
            raw_data[modulator] = raw_data[modulator].astype(int)
            raw_data = raw_data[(raw_data[variants[0]].isin([v1_case, v1_control])) & (raw_data[variants[1]].isin([v2_case, v2_control]))]
            # Plot for Binary traits
            if raw_data[trait].nunique() == 2:
                fig = binary_IATE_plot(raw_data, bqtl, modulator, trait)   
            # Plot for Continuous traits
            else:
                fig = continuous_IATE_plot(raw_data, bqtl, modulator, trait)

    elif param_type == "ATE":
        if "_&_" in treatments_string:
            st.markdown("Modulation plot can only be displayed when the Treatment contains only one variable.")
            return
        else:
            raw_data = pd.read_csv(
                raw_data_file(st.session_state.nextflow_rundir), 
                usecols = [trait, bqtl],
            ).dropna()
            raw_data[bqtl] = raw_data[bqtl].astype(int)
            case = float(selected["CASE"])
            control = float(selected["CONTROL"])
            raw_data = raw_data[(raw_data[bqtl].isin([case, control]))]
            # Plot for Binary traits
            if raw_data[trait].nunique() == 2:
                fig = binary_ATE_plot(raw_data, bqtl, trait)
            # Plot for Continuous traits
            else:
                fig = continuous_ATE_plot(raw_data, bqtl, trait)
        
    st.plotly_chart(fig, use_container_width=True)"""

def location_from_str(location_str):
    chr, start_end = location_str.split(":")
    start, end = start_end.split("-")
    return chr, int(start), int(end)

@st.cache_data
def SNPinfo(rsid, bqtls_data):
    response = http_variant_info(rsid)
    st.write(response)
    mapping_1 = response["mappings"][0]
    variant_row = bqtls_data[bqtls_data.ID == rsid].iloc[0]
    st.write(bqtls_data[bqtls_data.ID == rsid])
    chr, start, end = location_from_str(mapping_1["location"])
    if variant_row["REF.counts"] > variant_row["ALT.counts"]:
        binding_allele = variant_row.REF + " (" + str(variant_row["REF.counts"]) + ")"
        non_binding_allele = variant_row.ALT + " (" + str(variant_row["ALT.counts"]) + ")"
    else:
        non_binding_allele = variant_row.REF + " (" + str(variant_row["REF.counts"]) + ")"
        binding_allele = variant_row.ALT + " (" + str(variant_row["ALT.counts"]) + ")"
    
    basesnpinfo = {
        "Binding Allele (Counts)": binding_allele,
        "Non-Binding Allele (Counts)": non_binding_allele,
        "REF Allele": variant_row.REF,
        "ALT Allele": variant_row.ALT,
        "Chromosome": [chr],
        "Start": [start],
        "End": [end],
        "Strand": [mapping_1["strand"]],
        "Ensembl Alleles": [mapping_1["allele_string"]],
        "Ensembl Minor Allele": [response["minor_allele"]],
        "Ensembl MAF": [response["MAF"]],
        "Ensembl Ancestral Allele": [mapping_1["ancestral_allele"]],
    }
    st.header(f"Base Report")
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
def region_annotations(chr, v_start, v_end, distance, features):
    response = http_ensemble_annotations(
        chr, v_start, v_end, 
        distance=distance, 
        features=features
        )
    return pd.DataFrame(response)

@st.cache_data
def plot_motif_logo(json_response, pmin, pmax):
    bm_df = pd.DataFrame([json_response['elements'][str(i)] 
                          for i in range(1,json_response['length'])])
    logo = logomaker.Logo(bm_df)
    logo.highlight_position_range(pmin=pmin, pmax=pmax, color='gold', alpha=.5)
    logo.ax.set_title(
        f"Associated transcription factor complexes: "
        f"{','.join(json_response['associated_transcription_factor_complexes'])} "
        f"from {json_response['source']}"
    )
    st.pyplot(logo.fig)