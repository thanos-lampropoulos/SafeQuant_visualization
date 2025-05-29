## Title

SafeQuant_visualization

## Description

A Streamlit app that creates interactive volcano plots for SafeQuant protein quantitation data.

The Mass-Spec raw data from a label-free Proteomics experiment are first imported in Progenesis for Proteomics. The quantitation data from Progenesis are then imported in SafeQuant for statistical analysis.

The SafeQuant protein quantitation data can be upload in this Streamlit app, which will generate interactive volcano plots that depict the relative abundance of proteins between experimental conditions.


## How to Install and Run the Project

- All functions are written in Python 3.

- All dependencies are included in the .py file.

- The app has been deployed on [Streamlit Community Cloud](https://safequantvisualization.streamlit.app/). Alternatively, the .py file can be downloaded and executed using a local Streamlit installation.
  

# How to Use the Project

Use the provided PROTEIN.tsv file and upload it on the Streamlit app.

# Known issues

When generating your own data with SafeQuant, the experimental condition names should not contain any underscores.
