# TarViz
Visualization dashboard for TarGene results

## Development from Eddie

In order to forward and display the the app in the browser (workaround taken from the [jupyther notebook config](https://docs.anaconda.com/anaconda/user-guide/tasks/remote-jupyter-notebook/)):

1. Create a SSH tunnel from your local machine:

```bash
ssh -L 8080:localhost:8080 wild-kb
```
where wild-kb is the node where the app will be running and is described in the ssh config file.

2. Run the app:

```bash
streamlit run 0_🌶_Home.py --server.port=8080
```