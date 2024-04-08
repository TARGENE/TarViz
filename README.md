# TarViz
Visualization dashboard for TarGene results

## Requirements

Conda environment with poetry

## Development from Eddie

In order to forward and display the the app in the browser (workaround taken from the [jupyther notebook config](https://docs.anaconda.com/anaconda/user-guide/tasks/remote-jupyter-notebook/)):

1. Create a SSH tunnel from your local machine:

```bash
ssh -fNL 8501:localhost:8501 <Eddie node>
```

or

```
ssh -fN <UUN>@<Eddie_node> -J <UUN>@eddie.ecdf.ed.ac.uk -L 8501:localhost:8501
```
where <Eddie_node> is the node where the app will be running and is described in the ssh config file.

2. Run the app:

```bash
poetry run streamlit run 0_🌶_Home.py --server.port=8501 NEXTFLOW_RUNDIR <folder_to_GTEX/bulk-qtl_v8_multi-tissue-qtl_GTEx_Analysis_v8.hdf5>
```

where `NEXTFLOW_RUNDIR` is the TarGene pipeline's run directory.
