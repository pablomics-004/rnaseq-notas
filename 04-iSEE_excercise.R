## SummaryExperiment extension for scRNA-seq data
library("iSEE")

## Downloading data from spatialLIBD
sce_layer <- spatialLIBD::fetch_data("sce_layer")
iSEE::iSEE(sce_layer) # Exploring data with iSEE
