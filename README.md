# TarViz
Visualization dashboard for TarGene results

## Development from Eddie

In order to forward and display the the app in the browser (workaround taken from the [jupyther notebook config](https://docs.anaconda.com/anaconda/user-guide/tasks/remote-jupyter-notebook/)):

1. Create a SSH tunnel from your local machine:

```bash
ssh -fNL 8181:localhost:8181 wild-kb
```

or

```
ssh -fN s2042526@node2c17 -J s2042526@eddie.ecdf.ed.ac.uk -L 9999:localhost:9999 
```
where wild-kb is the node where the app will be running and is described in the ssh config file.

2. Run the app on Eddie
   Install poetry in conda environemnt on Eddie. Navigate to TarViz directory and run following commands
   
```bash
poetry init
poetry install
poetry run streamlit run 0_🌶_Home.py --server.port=8181 NEXTFLOW_RUNDIR
```
where `NEXTFLOW_RUNDIR` is the TarGene pipeline's run directory.

Server ports for local machine and Eddie must be the same. Run ```localhost:8181``` in browser to access dashboard. 
