#' create_probMatrix
#'
#' create_probMatrix creates a variable that contains the transition probabilities ("initial", "forward", "backward", "intermittent", "observed")
#'
#' @param input_data A dataset in a format similar to `initial_data`. This dataset must contain the variables "state_from", which is the status at the beginning of the transition (say smoker in 2010),
#' "state_to", which is the status at the end of the transition (say ex-smoker in 2011) and "tran_Year", which is an integer variable that is equal
#' to the number of transitions. "tran_Year" == 1 means that the transition occurs from 2010 to 2011, "tran_Year" == 2, from 2011 to 2012, up to
#' the total number of transitions
#' @param patient_id A character variable that specifies the column name with the unique Id of the patient
#'
#' @importFrom stats aggregate
#' @return a data frame containing the column "prob_matrix"
#' @export
#'
#' @examples
#' create_probMatrix(initial_data, patient_id = "patient_id")



create_probMatrix <- function(input_data, patient_id){

  input_data <- as.data.frame(input_data)
  stateStatuses <- data.frame(patient_id = input_data[, patient_id], tran_Year = input_data$tran_Year, state_to = as.numeric(!is.na(input_data$state_to)), state_from = as.numeric(!is.na(input_data$state_from)))
  stateStatuses$sumFromTo <- stateStatuses$state_from + stateStatuses$state_to
  numNonMissing <- aggregate(sumFromTo ~ patient_id, sum, data = stateStatuses)
  subj_initial <- numNonMissing$patient_id[which(numNonMissing$sumFromTo == 0)]
  initial_patients = input_data[which(input_data[, patient_id] %in% subj_initial),]
  if (nrow(initial_patients)>0){
    initial_patients$prob_matrix <- "forward"
    initial_patients$prob_matrix[which(initial_patients$prob_matrix=="forward" &
                                          initial_patients$tran_Year == 1)] <- "initial"
  }
  sub_back_int <- numNonMissing$patient_id[which(numNonMissing$sumFromTo >= 1)]
  back_interm_patients = input_data[which(input_data[, patient_id] %in% sub_back_int),]
  back_interm_patients$prob_matrix <- "nothing"

  idx.subj=unique(back_interm_patients[, patient_id])
  data <- data.frame()
  for (i in 1:length(idx.subj)){
    input_data <- back_interm_patients[back_interm_patients[, patient_id]==idx.subj[i],]
    rows <- nrow(input_data)

    idx.stateto = which(!is.na(input_data$state_to==T));idx.stateto
    idx.statefrom = which(!is.na(input_data$state_from==T));idx.statefrom
    id.all <- unique(sort(c(idx.stateto, idx.statefrom)));id.all

    if (length(id.all) > 1){
      for (iROW in 1:length(id.all)){
        if (iROW == length(id.all)) break
        input_data$prob_matrix[id.all[iROW]:id.all[iROW+1]] = "intermittent"
      }
    }

    id.first <- id.all[1]
    if (!is.na(input_data$state_to[id.first]) == T){
      input_data$prob_matrix[id.first:1] = "backward"
    }
    id.last <- id.all[length(id.all)]
    if (!is.na(input_data$state_from[id.last]) == T){
      input_data$prob_matrix[id.last:nrow(input_data)] = "forward"
    }

    data <- rbind(data,input_data)
  }

  final_data <- rbind(initial_patients, data)
  final_data$prob_matrix <- ifelse(is.na(final_data$state_from)==F & is.na(final_data$state_to)==F, "observed", final_data$prob_matrix)
  final_data <- final_data[order(final_data[, patient_id], final_data$tran_Year),]

  return(final_data)

}
