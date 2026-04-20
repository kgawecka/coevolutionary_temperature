# Geographic mosaic of coevolutionary temperature: from hotspots to coldspots in mutualistic communities

This repository contains the data and code used in the manuscript [Gawecka, K.A, Pedraza, F., Andreazzi, C.S., Bascompte, J. (2025) "Quantifying the geographic mosaic of coevolutionary temperature: from coldspots to hotspots"](https://www.biorxiv.org/content/10.1101/2025.11.27.691007v1.abstract).

## `Code`
All code was created in R version 4.5.1.

`preprocessing.R` - compute local network metrics

`coevolution_simulation.R` - simulate mutualistic coevolution on local networks

`analysis.R` - analyse coevolution model output and plot manuscript figures

## `Data`
Local interaction networks and patch data were kindly provided by Frank Jauker, Birgit Jauker and Ingo Grass ([Grass et al., 2018](https://www.nature.com/articles/s41559-018-0631-2); [Jauker et al., 2019](https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecy.2569)).

`data_networks.csv` - local network metrics (output of `preprocessing.R`)

`data_patches.csv` - habitat patch data

`data_species.csv` - local species degree (output of `preprocessing.R`)

`Minc_XX.csv` - local interaction network incidence matrix (XX - network number)

## `Output`

`df_RS_XX_mXX_alphaXX.csv` - coevolutionary model result (output of `coevolution_simulation.R`; XX - network number / m-value / alpha-value)
