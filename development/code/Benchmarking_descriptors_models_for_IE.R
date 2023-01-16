library(tidyverse)
library(caret)
library(caTools)
library(Metrics)
library(plotly)
library(cowplot)
setwd("C:/Users/HelenSepman/OneDrive - Kruvelab/Documents/GitHub/MS2Quant/development")
source("code/functions.R")
library("FactoMineR")
library("factoextra")

# ---- theme ----
extrafont::loadfonts(device = "win")

font <- extrafont::choose_font("Quicksand")
fontsize <- 14
basecolor <- "#515251" 
highlighter_color = "#7CB368" #"#728D68"
datapoints_color1 <- "#959D95"
datapoints_color2 <- "#728D68"

my_theme <-   theme(
  plot.background = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(size = 1,
                           color = basecolor),
  # plot.title = element_text(color = basecolor,
  #                           size = 14,
  #                           face = "bold"),
  text = element_text(family = font,
                      size = fontsize,
                      color = basecolor),
  legend.key = element_blank(),
  strip.background = element_blank(),
  strip.text = element_blank(),
  # strip.text = element_text(family = font,
  #                           size = fontsize,
  #                           color = basecolor),
  legend.position = "none",
  legend.title = element_blank(), 
  legend.text = element_text(family = font,
                             size = fontsize,
                             color = basecolor),
  axis.text = element_text(family = font,
                           size = fontsize,
                           color = basecolor),
  axis.ticks = element_blank(),
  aspect.ratio = 1,
  axis.title.x = element_text(hjust = c(1), vjust = c(0)),
  axis.title.y = element_text(hjust = c(1), vjust = c(1))
)

# ---- Read in data, create split for train and test set -----

data <- read.csv("C:/Users/HelenSepman/OneDrive - Kruvelab/Documents/GitHub/MS2Quant/development/data/data_unified_IEs_all_filtered_20221121.csv")

inchi <- data %>%
  select(inchi) %>%
  unique()

set.seed(123)
split_inchi <- sample.split(inchi$inchi, SplitRatio = 0.8)

inchi <- inchi %>%
  mutate(split = split_inchi)

data <- data %>%
  left_join(inchi)

# Distribution of data
ggplot(data = data) + 
  geom_histogram(mapping = aes(x=as.numeric(unified_IEs)), bins = 45) +
  facet_wrap(~split, scales="free_y") +
  xlab(substitute(paste("log", italic("IE"))["unified"]))  +
  my_theme

# -----------------------------------------------------
# Calculating (reading in) descriptors, cleaning them 
# -----------------------------------------------------

setwd("C:/Users/HelenSepman/OneDrive - Kruvelab/Documents/GitHub/MS2Quant/development/code")

# ---- Cleaning SMILES ----
data_SMILES <- data %>%
  select(SMILES) %>%
  unique()
  
# ----  PaDEL descriptors ---- 
descr_padel <- PaDEL_original(data_SMILES)

# ----  Mordred descriptors ----

descr_mordred <- Mordred_descriptors(data_SMILES)

# ----  Morgan-2 descriptors ----

descr_morgan <- Morgan2_descriptors(data_SMILES)

# ----  MAP4 ----
#calculated on a Linux-based VM, read in here:
descr_map4 <- read_delim("MAP4fingerprints.csv", delim = ",", col_names = FALSE)
smiles_map4 <- read_delim("MAP4fingerprints_smiles.csv", delim = ",", col_names = FALSE)

colnames(smiles_map4) <- c("SMILES")

descr_map4 <- smiles_map4 %>%
  bind_cols(descr_map4)

descr_map4 <- data_SMILES %>% 
  left_join(descr_map4)


# ---- rcdk ----

descr_rcdk <- rcdk_fingerprints(data_SMILES)

# ---- SIRIUS structural fingerprints ----

SIRIUS_descriptors <- Fingerprint_calc(data_SMILES)

# NB! filter out the FPs that SIRIUS actually does not calculate!!
sirius_true_FP_pos <- read_delim("csi_fingerid.tsv",
                                 delim = "\t")
sirius_true_FP_neg <- read_delim("csi_fingerid_neg.tsv",
                                 delim = "\t")


col_names_calcSirius_pos <- as.vector(paste("Un", sirius_true_FP_pos$absoluteIndex, sep = ""))
col_names_calcSirius_neg <- as.vector(paste("Un", sirius_true_FP_neg$absoluteIndex, sep = ""))
col_names_calcSmiles <- colnames(SIRIUS_descriptors)
intersection <- intersect(col_names_calcSirius_pos, col_names_calcSirius_neg)
intersection <- intersect(intersection, col_names_calcSmiles)
intersection <- append("SMILES", intersection)

SIRIUS_descriptors <- SIRIUS_descriptors[names(SIRIUS_descriptors)[names(SIRIUS_descriptors) %in% intersection] ]


# ---- Cleaning descriptors ----

descr <- SIRIUS_descriptors

descr <- cleaning_descriptors(descr)

# -----------------------------------------------------
# Split data to training and test set, building models
# -----------------------------------------------------

data_clean <- data %>%
  select(unified_IEs, inchi, split, SMILES, viscosity, surface_tension, polarity_index, pH_aq, name, Lab) %>% 
  left_join(descr) %>%
  unique() %>%
  na.exclude()

training_set = data_clean %>%
  filter(split) %>%
  select(-split)

test_set = data_clean %>%
  filter(!split) %>%
  select(-split)

# ---- Cross-validation/hyperparameter tuning ----
set.seed(123)
folds = 5
fitControlMethod = "boot"
method = "xgbTree"             
folds = groupKFold(training_set$inchi, k = folds) 
fitControl <- trainControl(method = fitControlMethod, index = folds)

# tested models
# 1) xgbTree
# 2) xgbLinear
# 3) xgbDART


# ---- Training a model ----
set.seed(123)
model <- train(unified_IEs ~ ., 
               data = training_set %>% select(-inchi, -SMILES, -name, -Lab),
               method = method,
               trControl = fitControl)

# --------


# Predict logIE values for training set; RMSE of training set
training_set = training_set %>%
  mutate(unified_IEs_pred = predict(model, newdata = training_set))

rmse_training <- rmse(training_set$unified_IEs, training_set$unified_IEs_pred)

# Predict logIE values for test set, RMSE of test set
test_set <- test_set %>%
  mutate(unified_IEs_pred = predict(model, newdata = test_set))

rmse_test <- rmse(test_set$unified_IEs, test_set$unified_IEs_pred)


# Plot predicted vs experimental values for training and test set

data_forPlotting <- training_set %>%
  mutate(set = "training") %>%
  bind_rows(test_set %>%
              mutate(set = "test")) %>%
  mutate(set = factor(set,
                      ordered = TRUE,
                      levels = c( "training", "test")))


plot = ggplot(data = data_forPlotting) + 
  geom_point(
             mapping = aes(x = unified_IEs,
                           y = unified_IEs_pred,
                           color = set,
                           text = name),
             size = 2.5,
             alpha = 0.5) +

  scale_color_manual(values=c("#515251", "#7CB368"))+

  geom_abline(intercept = 0, slope = 1) +
  geom_abline(intercept = 1, slope = 1) +
  geom_abline(intercept = -1, slope = 1) +
  theme(aspect.ratio = 1) + 
  xlim(c(-2.5, 6.5)) +
  ylim(c(-2.5, 6.5)) +
  #annotation_logticks() +
  #theme_plot +
  xlab(substitute(paste("log", italic("IE"))["unified"]))  +
  ylab(substitute(paste("log", italic("IE"))["predicted"])) +
  my_theme #+
  #facet_wrap(~ set, scales = "fixed")
  
plot

data_all <- list("model" = model,
                 "plot" = plot,
                 "data" = data_clean,
                 "data_forPlotting" = data_forPlotting,
                 "rmse_training" = rmse_training,
                 "rmse_test" = rmse_test,
                 "date" = Sys.Date())


# ---- Save data out as list ----

saveRDS(data_all, file="model_MAP4_xgbDART_test.RData")




# ---- Interpretation of results: statistics and prediction error calculation ----

# Q2 for the training set
mean(model$resample$Rsquared)

# R2 for the test set
r = cor(data_forPlotting$unified_IEs, data_forPlotting$unified_IEs_pred)
r_squared = r^2

# mean, median, geom. mean
test_set_pred_error <- data_forPlotting %>%
  filter(set == "test") %>% 
  mutate(pred_error = case_when(unified_IEs > unified_IEs_pred ~ (10^unified_IEs)/(10^unified_IEs_pred),
                                TRUE ~ (10^unified_IEs_pred)/(10^unified_IEs))) %>%
  select(pred_error, unified_IEs, unified_IEs_pred, everything())

# mean pred error
mean(test_set_pred_error$pred_error)

# geometric mean
exp(mean(log(test_set_pred_error$pred_error)))

# median pred error
median(test_set_pred_error$pred_error)

summary(lm(unified_IEs_pred ~ unified_IEs, data = test_set_pred_error))
# statistically significatn - define alpha, in M&M





# ---- training for best fit -----
training_set = data_clean %>%
  select(-split)

grid_optimized <- expand.grid(
  nrounds = model$bestTune$nrounds,
  max_depth = model$bestTune$max_depth,
  eta = model$bestTune$eta,
  gamma = model$bestTune$gamma,
  colsample_bytree = model$bestTune$colsample_bytree,
  min_child_weight = model$bestTune$min_child_weight,
  subsample = model$bestTune$subsample
)

train_control <- caret::trainControl(
  method = "none"
)

set.seed(123)
model2 <- train(unified_IEs ~ .,
                method = method,
                data = training_set %>% select(-inchi, -SMILES, -name, -Lab),
                trControl = train_control,
                tuneGrid = grid_optimized)

# Predict logIE values for training set, RMSE of training set
training_set = training_set %>%
  mutate(unified_IEs_pred = predict(model2, newdata = training_set))

rmse_training <- rmse(training_set$unified_IEs, training_set$unified_IEs_pred)


data_all <- list("model" = model2,
                 "data" = data_clean,
                 "rmse_training" = rmse_training,
                 "date" = Sys.Date())

# ---- Save data out as list ----
# saveRDS(data_all, file="model_PaDEL_xgbTree_allData.RData")



