# MS2Quant

MS2Quant package has been tested with R version 4.1.1 and with SIRIUS versions 4.9.15 to 5.6.2.

To install the package:
```
devtools::install_github("kruvelab/MS2Quant"
                         ,ref="main"
                         ,INSTALL_opts="--no-multiarch")
```


Test code using example data from the package:
```
library(MS2Quant)

path_dataframe_calibrants_suspects <- system.file("example_data", "quantification_example.csv", package = "MS2Quant")
path_eluent_file <- system.file("example_data", "eluent.csv", package = "MS2Quant")
path_suspects_sirius_project_folder <- system.file("example_data", "SIRIUS_results", package = "MS2Quant")

MS2Quant_quantification_results <- MS2Quant_quantify(path_dataframe_calibrants_suspects,
                                                     path_eluent_file ,
                                                     organic_modifier = "MeCN",
                                                     pH_aq = 2.7,
                                                     path_suspects_sirius_project_folder)

# Separate calibration plots for each calibrant
MS2Quant_quantification_results$calibrants_separate_plots

# Calibration plot between experimental logRF and logIE
MS2Quant_quantification_results$logIE_logRF_calibration_plot

# Summary of the linear model between logRF and logIE
MS2Quant_quantification_results$calibration_linear_model_summary

# All suspect concentrations
MS2Quant_quantification_results$suspects_concentrations

# Date when the quantification was done
MS2Quant_quantification_results$date


```

