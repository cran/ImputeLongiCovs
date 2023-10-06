#' Initial Data for imputing categorical covariates
#'
#' A dataset containing longitudinal data. The outcome of interest is the smoking status with three states (smoker, exsmoker, neversmoker),
#' which are represented via transitions.
#'
#' \itemize{
#'   \item {patient_id}: {Unique identifier for each patient}
#'   \item {tran_Year}: {numeric, starting from 1 up to the number of transitions}
#'   \item {transition_year}: {text explanation of the transition}
#'   \item {state_from}: {the state at the beginning of a transition}
#'   \item {state_to}: {the state at the end of a transition}
#'   \item {cardio_state_from}: {cardiovascular disease at the beginning of the transition, binary, if 1 == Yes, else No}
#'   \item {cardio_state_to}: {cardiovascular disease at the end of the transition, binary, if 1 == Yes, else No}
#'   \item {flu_vaccination_state_from}: {flu vaccination at the end of the transition, binary, if 1 == Yes, else No}
#'   \item {flu_vaccination_state_to}: {flu vaccination disease at the end of the transition, binary, if 1 == Yes, else No}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name initial_data
#' @usage data(initial_data)
#' @format A data frame with 2000 rows and 9 variables
NULL

