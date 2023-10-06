#' impute_categorical_covariates
#'
#' impute_categorical_covariates imputes longitudinal categorical covariates through a joint model that accommodates
#' initial, forward, backward, and intermittent transitions.
#'
#' It encloses three different functions. The `initial_forward_function` imputes the longitudinal categorical covariate of interest based on whether in that transition the `prob_matrix` of a patient was `initial` or `forward`.
#' The `imputeIntermittent` imputes the longitudinal categorical covariate for the intermittent transition and the `backward_function` imputes the longitudinal categorical covariate for the backward transition.
#'
#' @param input_data A dataset in a format similar to `analyses_data`. This dataset must contain the variables "state_from", which is the status at the beginning of the transition (say smoker in 2010),
#' "state_to", which is the status at the end of the transition (say ex-smoker in 2011) and "tran_Year", which is an integer variable that is equal
#' to the number of transitions. "tran_Year" == 1 means that the transition occurs from 2010 to 2011, "tran_Year" == 2, from 2011 to 2012, up to
#' the total number of transitions Also, it must contain "prob_matrix" which captures all the transitions ("initial", "forward", "backward", "intermittent", "observed") that was
#' calculated with the `create_probMatrix` function
#' @param patient_id A character variable that specifies the column name with the unique Id of the patient
#' @param number_of_transitions The number of transitions needed. For example for years 2010, 2011 and 2012 there exist 2 transitions.
#' @param covariates_initial The covariates to be used in the initial model
#' @param covariates_transition The covariates to be used in the transition model
#' @param missing_variable_levels The levels of the missing categorical outcome (e.g. "smoker", "ex-smoker", "never-smoker")
#' @param startingyear If the starting year per patient has no missing values, specify it
#' @param without_trans_prob This statement is useful when there are very high proportions of missing data and our initial and transition model cannot converge.
#' It provides the user with two options. One, to "notImpute", namely to return NA and two, to "ImputeEqualProbabilities", i.e., the user
#' can sample with equal probabilities.
#' @param m Numeric, the number of imputed datasets
#'
#' @importFrom nnet multinom
#' @importFrom stats complete.cases as.formula predict

#' @return a list of m data frames with no missing values in the categorical outcome
#' @export
#'
#' @references ()
#'
#' @examples
#' impute_categorical_covariates(analyses_data,
#' patient_id = "patient_id",
#' number_of_transitions = 2,
#' covariates_initial = c("cardio_state_from", "flu_vaccination_state_from"),
#' covariates_transition = c("cardio_state_to", "flu_vaccination_state_to"),
#' missing_variable_levels = c("never-smoker", "smoker", "ex-smoker"),
#' startingyear = NULL,
#' without_trans_prob = "notImpute",
#' m = 2)


# Then I can do multiple imputations
impute_categorical_covariates <- function(input_data,
                                          patient_id,
                                          number_of_transitions,
                                          covariates_initial = NULL,
                                          covariates_transition = NULL,
                                          missing_variable_levels,
                                          startingyear = NULL,
                                          without_trans_prob ,
                                          m = 1){

  imputed_datasets <- list()

  for (iRep in 1:m) {
    imputed_datasets[[iRep]] <- impute_categorical_covariates_once(input_data = input_data,
                                                                   patient_id = patient_id,
                                                                   number_of_transitions,
                                                                   covariates_initial = covariates_initial,
                                                                   covariates_transition = covariates_transition,
                                                                   missing_variable_levels = missing_variable_levels,
                                                                   startingyear = NULL,
                                                                   without_trans_prob = without_trans_prob)
    imputed_datasets[[iRep]]$input_data$.imp <- iRep
  }

  return(imputed_datasets)
}

