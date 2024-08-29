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

2. Run the app:

```bash
streamlit run 0_🌶_Home.py --server.port=6666 NEXTFLOW_RUNDIR
```

where `NEXTFLOW_RUNDIR` is the TarGene pipeline's run directory.


## Information about running TarViz from a local machine

Required informaton (examples found in example folder):
1. A CSV with information about results, including p-values, treatments, and outcomes
2. CSV with information about bQTLs from baal-nf pipeline including CHR, POS, RSID, REF, ALT, REF.counts, ALT.counts
3. Nextflow config file that points to results and bQTL information.
4. Data from GTeX - this can be found on the UKBB-53116 folder on datastore.

  ### Steps to run TarViZ locally

1. Create poetry environment **conda env create -n poetry poetry**
2. Activate poetry environment **conda activate poetry**
3. Initialise poetry with **poetry init**
4. Save a folder with above information somewhere locally
5. Navigate to *tarviz/cwd/tarviz_test* folder which should have a *pyproject.toml* in it and **run poetry install**
6. Ensure all packages have installed correctly (if on Macbook you may have to do a brew install hdf5 to make tables install but any other issues )
7. Make sure connection to Datastore is active (For Macbook users: Open Finder, CMD + K, if you nothing there type in *smb://cmvm.datastore.ed.ac.uk/igmm* and follow instructions, if smb://cmvm.datastore.ed.ac.uk/igmm active then proceed)
8. **Open Terminal** and navigate to the correct folder and run the **command poetry run streamlit run 0_Home.py ../data/tarviz_folder/ /Volumes/igmm/UK-BioBank-53116/other/gtex_data/chrom/**
