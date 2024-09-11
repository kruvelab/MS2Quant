# MS2Quant

**MS2Quant** is a machine learning model that enables prediction of concentration 
from fragmentation (MS<sup>2</sup>) spectra of detected, but unidentified chemicals. 
MS2Quant is an *xgbTree* algorithm-based regression model developed using ionization efficiency data for
1191 unique chemicals that span 8 orders of magnitude. For more information about the data used in model training 
see the summary file under [MS2Quant/development/data](https://github.com/kruvelab/MS2Quant/tree/main/development/data).

As for now, MS2Quant can predict ionization efficiencies in positive electrospray ionization mode
for [M+H]+ and [M]+ ions, and the negative mode model is currently being developed.
MS2Quant package has been tested with R version 4.1.1 and 4.2.2, as well as with SIRIUS versions 4.9.15 to 5.6.2.

For more information about development and validation of MS2Quant, please see the [paper](https://doi.org/10.1021/acs.analchem.3c01744).

## 1. Using MS2Quant to quantify unidentified substances

Following steps need to be taken to use MS2Quant for quantification of detected unknown chemicals:

+ 1.1. Install the MS2Quant R-package
+ 1.2. Measure a set of calibrants together with the sample
+ 1.3. Organize the input information into a table
+ 1.4. Run SIRIUS+CSI:FingerID calculations
+ 1.5. Use MS2Quant to quantify

![MS2Quant_github_workflow](https://github.com/kruvelab/MS2Quant/assets/48623628/1d32b3bc-8cd0-4ecf-99c1-d2630f6f5af6)
**Figure 1.** MS2Quant workflow


### 1.1. Install the MS2Quant R-package

To install the package, run the following command in R:

```
devtools::install_github("kruvelab/MS2Quant",
                         ref="main",
                         INSTALL_opts="--no-multiarch")
```
To test if the installation was successful, you can run the test code (see testing in Chapter 4.)

### 1.2. Measure a set of calibrants together with the sample

To convert a predicted ionization efficiency (IE) value to an instrument and measurement-specific response 
factor (RF), calibration of the model is performed by measuring calibrants during the same experimental run with suspects.
Response factors of calibrants are calculated by MS2Quant using spiked concentration and integrated signal areas, and the 
linearity is checked (see *1.2.1. Linearity check*). 
To construct a calibration graph between predicted ionization efficiencies and experimental response factors, 
the ionization efficiency of calibrants is predicted using their SMILES notation, and the calibration
graph between logIE and logRF can be visualized (see testing in Chapter 4.).

#### 1.2.1. Linearity check
Weighted linear regression is performed for calibration curves (spiked concentration vs integrated signal area). 
The goodness of fit evaluated based on relative residuals, and if the largest relative residual 
(based on absolute value) is >20%, highest concentration datapoint is removed. If this is not achieved, minimum 5 points are 
left in the dataset, however, a warning message will be sent if the largest relative residual is larger than 20%. We suggest looking
into the calibration graphs visually to confirm the correctness of provided concentrations and signal areas 
(separate calibration plots for each calibrant can be visualized after running *MS2Quant_quantify()* function, see testing in Chapter 4.)

#### 1.2.2. Choosing calibrants
The calibrants are used to find a linear regression between predicted ionization effciencies and measurement
specific response factors. Therefore, it is important that the ionization efficiencies of these chemicals would cover 
multiple orders of magnitude in IE scale. Minimum 5 calibrants is required to use MS2Quant but we strongly suggest to use 
at least 10 calibrants. For choosing calibrants that would cover wide range of ionization efficiencies, you can look into experimental IE data used for
model training at [MS2Quant/development/data](https://github.com/kruvelab/MS2Quant/tree/main/development/data).

### 1.3. Organize the input information into a table
To use the function *MS2Quant_quantify()*, all data about calibrants and suspects needs to be organized into a table
with following information:
+ *identifier* - this information is especially relevant for unidentified features as the identifier needs to align with the identifier in 
SIRIUS+CSI:FingerID calculations folder name (FIG). For calibrants, this information is used for grouping; therefore, the name of the chemical can be added, for example.
+ *SMILES* - SMILES notation needs to be added to calibrants, as for these chemicals ionization efficiency is predicted from the structure. 
+ *retention_time* - retention time information is needed for all features (both calibrants and suspects subject to quantification). This information is used to
calculate eluent descriptors such as organic modifier percentage, aqueous pH, polarity index, surface tension, and viscosity.
+ *area* - integrated signal area; needs to be added for all LC/MS features. For calibrants, all signal areas corresponding to different spiked concentration levels should be added as a new line.
+ *conc_M* - molar concentration; needs to be added to calibrants in order to fit calibration curves and calculate RFs. **This needs to be left empty for chemicals subject to quantification!**

NB! Do not change the writing of the column names!
The example data table can be found at [inst/example_data/quantification_example.csv](https://github.com/kruvelab/MS2Quant/blob/main/inst/example_data/quantification_example.csv)


### 1.4. Run SIRIUS+CSI:FingerID
SIRIUS calculations can be run either using graphical user interface or command line interface.
Information obtained from SIRIUS predictions will be joined to the unidentified chemicals information based on the **identifier**. 
That means that the identifier provided in the table for unidentified LC/MS features needs to align with identifier found in SIRIUS folder name, see Figure 2.

**NB! Do not use underscores (_) in your identifier!**

![MS2Quant_github_sirius](https://github.com/kruvelab/MS2Quant/assets/48623628/89874378-64bc-499e-aa32-43c9f69447b1)
**Figure 2.** *identifier* of unidentified suspect will be used to join predicted fingerprints from SIRIUS calculations.


Longer explanation on how to use SIRIUS for predictions can be found from [Official online documentation for the SIRIUS MS/MS software github page](https://boecker-lab.github.io/docs.sirius.github.io/)
and some examples can be found on [MS2Tox github page](https://github.com/kruvelab/MS2Tox/tree/main#fingerprint-calculation-with-siriuscsifingerid).

### 1.5. Use MS2Quant to quantify
To quantify the unidentified chemicals, *MS2Quant_quantify()* function can be used. Following inputs are neede for the function:
+ *calibrants_suspects* - (path to) data table including information about calibrants and suspects subject to quantification (see Chapter 1.3.)
+ *eluent* - (path to) the gradient program information (see example at [inst/example_data/eluent.csv](https://github.com/kruvelab/MS2Quant/blob/main/inst/example_data/eluent.csv))
+ *organic_modifier* - specify which organic modifier was used (either "MeCN" or "MeOH").
+ *pH_aq* - specify the pH of the water phase.
+ *fingerprints* - specify the path to your SIRIUS calculations results folder.


## 2. Using MS2Quant to quantify (tentatively) identified substances

MS2Quant can also be used to quantify identifier or tentatively identified structures. For this, the workflow is similar to 
Chapter 1, with a modification in creating the input table (Chapter 1.3.). For tentatively identified chemicals, SMILES notation needs to be added 
(however, *conc_M* should still be left empty to consider this chemical subject to quantification).

![MS2Quant_github_tentativelyidentified](https://github.com/kruvelab/MS2Quant/assets/48623628/018c8a8d-443f-4274-ab6b-099c761d3b64)
**Figure 3.** If SMILES notation is provided for tentatively identified suspect, quantification will be done based on structure.

**NB!** Unidentified and tentatively identified chemicals subject to quantification can be used at the same time; however, when SMILES is provided,
quantification based on structure will be prioritized.


## 3. Using MS2Quant to predict ionization efficiency (without quantification) 

Quantification cannot be performed if calibration chemicals have not been run together with the suspects. Still, ionization efficiency can be predicted for any chemical either based on structure (SMILES) or
predicted fingerprints.

In order to predict ionization efficiency values to different chemicals under the same eluent conditions, following code can be run:
```
chemicals_SMILES = tibble(SMILES = c("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", 
                                     "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"))

IE_pred = MS2Quant_predict_IE(chemicals_SMILES)
data = IE_pred$chemicals_predicted_IEs
```
By default, 80/20 MeCN/H2O with aqueous pH of 2.7 is used. Alternative eluent compositions can be specified as function input: 
```
IE_pred = MS2Quant_predict_IE(chemicals_for_IE_prediction = chemicals_SMILES,
                              organic_modifier = "MeOH",
                              organic_percentage = 50,
                              pH_aq = 7)
```
Ionization efficiency can also be predicted for unidentified chemicals. In this case, similarly to *MS2Quant_quantify()* function, data table with an identifier and retention time needs to be provided together
with the eluent file and SIRIUS calculations results folder:
```
IE_pred = MS2Quant_predict_IE(chemicals_for_IE_prediction = path_dataframe_calibrants_suspects,
                              eluent = path_eluent_file,
                              fingerprints = path_suspects_sirius_project_folder,
                              organic_modifier = "MeOH",
                              pH_aq = 7)
```
If retention time is missing from the data table, default eluent conditions will be used. 

## 4. Test code using example data from the package

To confirm that the installation of **MS2Quant** package has been successful, the following code can be run.
This code uses test data from the package.

```
library(MS2Quant)

# positive mode MS2Quant
path_dataframe_calibrants_suspects <- system.file("example_data", "quantification_example.csv", package = "MS2Quant")
path_eluent_file <- system.file("example_data", "eluent.csv", package = "MS2Quant")
path_suspects_sirius_project_folder <- system.file("example_data", "SIRIUS_results_pos", package = "MS2Quant")
ionization = "esi_pos"

# negative mode MS2Quant
path_dataframe_calibrants_suspects <- system.file("example_data", "quantification_example_neg.csv", package = "MS2Quant")
path_eluent_file <- system.file("example_data", "eluent.csv", package = "MS2Quant")
path_suspects_sirius_project_folder <- system.file("example_data", "SIRIUS_results_neg", package = "MS2Quant")
ionization = "esi_neg"

MS2Quant_quantification_results <- MS2Quant_quantify(path_dataframe_calibrants_suspects,
                                                     path_eluent_file ,
                                                     organic_modifier = "MeCN",
                                                     pH_aq = 2.7,
						     NH4 = 0,
						     model = "MS2Quant",
						     ionization,
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

Complementary models on the same datasets have been developed based on PaDEL descriptors (model = "PaDEL").
Additionally, [negative mode PFAS quantification model by Lauria et al.](https://doi.org/10.1021/acs.est.3c07220) can be accessed (model = "model_PFAS").

