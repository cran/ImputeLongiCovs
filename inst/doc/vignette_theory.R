## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("ImputeLongiCovs")

## -----------------------------------------------------------------------------
library(ImputeLongiCovs)

## -----------------------------------------------------------------------------
data(input_data = initial_data, package = "ImputeLongiCovs")
head(initial_data, 7)

## -----------------------------------------------------------------------------
data(analyses_data, package = "ImputeLongiCovs")
head(analyses_data, 7)

## -----------------------------------------------------------------------------
data(initial_data, package = "ImputeLongiCovs")
head(initial_data, 7)

## -----------------------------------------------------------------------------
initial_data_after_function <- create_probMatrix(initial_data, patient_id = "patient_id")
head(initial_data_after_function, 7)

## -----------------------------------------------------------------------------
initial_data_subset <- initial_data_after_function[which(initial_data_after_function$patient_id %in% c("patient_10", "patient_102", "patient_114", "patient_136")),]
initial_data_subset <- initial_data_subset[, c(1:5)]
initial_data_subset

## -----------------------------------------------------------------------------
data(analyses_data, package = "ImputeLongiCovs")
head(analyses_data, 7)

## ---- results='hide'----------------------------------------------------------
imputed_smoking_status  <- impute_categorical_covariates(input_data = analyses_data, 
                                                         patient_id = "patient_id", 
                                                         number_of_transitions = 2, 
                                                         covariates_initial = c("cardio_state_from", "flu_vaccination_state_from"),
                                                         covariates_transition = c("cardio_state_to", "flu_vaccination_state_to"),
                                                         missing_variable_levels = c("never-smoker", "smoker", "ex-smoker"), 
                                                         startingyear = NULL,
                                                         without_trans_prob = "notImpute",
                                                         m = 1)

## -----------------------------------------------------------------------------
imputed_smoking_status <- imputed_smoking_status[[1]]$input_data
imputed_smoking_status <- imputed_smoking_status[, c(1:5)]

head(imputed_smoking_status, 21)

## ----secondgraph, echo=FALSE, out.width = '70%', fig.align="center"-----------
knitr::include_graphics("Flowchart_long_imputation1.png")

