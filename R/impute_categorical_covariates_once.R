#' impute_categorical_covariates_once
#'
#' impute_categorical_covariates_once imputes categorical covariates through a joint model that accommodates
#' initial, forward, backward, and intermittent transitions
#'
#' It encloses three different functions. The first deals with the initial and forward transitions, the second
#' with the intemittent and the third with the backward.
#'
#' @param input_data a data frame that contains the variables "state_from", which is the status at the beginning of the transition (say smoker in 2010),
#' "state_to", which is the status at the end of the transition (say ex-smoker in 2011) and "tran_Year", which is an integer that is equal
#' to the number of transitions (say tran_Year = 1, means that the transition from 2010 to 2011, tran_Year = 2, from 2011 to 2012, up to
#' the total number of transitions). Also, "prob_matrix" contains all the transitions ("initial", "forward", "backward", "intermittent") that were
#' calculated with the extract_probMatrix function
#' @param patient_id A character variable that specifies the column name with the unique Id of the patient
#' @param number_of_transitions The number of transitions needed. For example for years 2010, 2011 and 2012 there are 2 transitions.
#' @param covariates_initial The covariates to be used in the initial model
#' @param covariates_transition The covariates to be used in the transition model
#' @param missing_variable_levels The levels of the missing categorical outcome (e.g. "smoker", "ex-smoker", "never-smoker")
#' @param startingyear If the starting year per patient has no missing values, specify it
#' @param without_trans_prob This statement is useful when there are very high proportions of missing data and our initial and transition model cannot converge.
#' It provides the user with two options. One, to "notImpute", namely to return NA and two, to "ImputeEqualProbabilities", i.e., the user
#' can sample with equal probabilities.
#'
#' @importFrom nnet multinom
#' @importFrom stats complete.cases as.formula predict

#' @return a data frame with no missing values in the categorical outcome
#' @noRd
#'
#'
#' @examples
#' impute_categorical_covariates_once(analysis_data,
#'  patient_id = "patient_id",
#'  number_of_transitions = 2,
#'  covariates_initial = c("cardio_state_from", "flu_vaccination_state_from"),
#'  covariates_transition = c("cardio_state_to", "Diabetes_state_to", "flu_vaccination_state_to"),
#'  missing_variable_levels = c("never-smoker", "ex-smoker", "smoker"),
#'  startingyear = startingyear,
#'  without_trans_prob = "notImpute")




impute_categorical_covariates_once <- function(input_data,
                                               patient_id,
                                               number_of_transitions,
                                               covariates_initial = NULL,
                                               covariates_transition = NULL,
                                               missing_variable_levels,
                                               startingyear = NULL,
                                               without_trans_prob = c("notImpute", "ImputeEqualProbabilities")){
  input_data <- as.data.frame(input_data)

  if(!without_trans_prob %in% c("notImpute", "ImputeEqualProbabilities")){
    stop("without_trans_prob should take notImpute or ImputeEqualProbabilities values")
  }

  if(!"state_from" %in% colnames(input_data)){
    stop("state_from column does not exist")
  }

  if(!"state_to" %in% colnames(input_data)){
    stop("state_to column does not exist")
  }

  if(!"tran_Year" %in% colnames(input_data)){
    stop("tran_Year column does not exist")
  }
  if(!is.numeric(input_data$tran_Year) == T){
    stop("tran_Year column should be integer")
  }
  input_data$state_from <- factor(as.character(input_data$state_from), levels = missing_variable_levels)
  input_data$state_to <- factor(as.character(input_data$state_to), levels = missing_variable_levels)

  unique_values_levels = unique(missing_variable_levels);unique_values_levels
  unique_values_state_from = unique(input_data$state_from);unique_values_state_from
  unique_values_state_to = unique(input_data$state_to);unique_values_state_to

  if(any(unique_values_levels %in% unique_values_state_from == F)){
    stop("All levels of missing_variable_levels vector should be present in the state_from column")
  }
  if(any(unique_values_levels %in% unique_values_state_to == F)){
    stop("All levels of missing_variable_levels vector should be present in the state_to column")
  }

  observed_smoke_transition <- input_data[complete.cases(input_data),]
  observed_smoke_initial <- input_data[complete.cases(input_data$state_from),]


  draw_sample <- function(x, missing_variable_levels, without_trans_prob){
    if (any(is.na(x))){
      if (without_trans_prob == "notImpute"){
        return(NA)
      } else {
        sample(missing_variable_levels, 1)
      }
    } else {
      sample(missing_variable_levels, 1, prob = x)
    }
  }


  if (is.null(covariates_initial)){
    formulainitial <- as.formula(paste("state_from ~ 1", sep = ""))
  } else {
    formulainitial <- as.formula(paste("state_from ~ ", paste(covariates_initial, collapse = "+"), sep = ""))
  }

  if (is.null(covariates_transition)){
    formulatransition <- as.formula(paste("state_to ~ state_from", sep = ""))
  } else {
    formulatransition <- as.formula(paste("state_to ~ state_from + ", paste(covariates_transition, collapse = "+"), sep = ""))
  }
  modelFittransition <- try(multinom(as.formula(formulatransition), data = observed_smoke_transition, trace = FALSE))
  if(is.character(modelFittransition)){
    stop("Error: The transition model could not be fitted")
  }

  modelFitinitial <- try(multinom(as.formula(formulainitial), data = observed_smoke_initial, trace = FALSE))
  if(is.character(modelFitinitial)){
    stop("Error: The initial model could not be fitted")
  }


  if (is.null(startingyear)& any(input_data$prob_matrix %in% "initial" == T)){


    #### 2.1 initial, forward  function ####
    initial_forward_function <- function(input_data = input_data,
                                         number_of_transitions = number_of_transitions,
                                         covariates_initial = covariates_initial,
                                         covariates_transition = covariates_transition,
                                         missing_variable_levels = missing_variable_levels){



      initialProbs <- predict(modelFitinitial, newdata = input_data[which(input_data$prob_matrix == "initial"), covariates_initial],"prob")
      check_data <- input_data[which(input_data$prob_matrix == "initial"), covariates_initial]

      if (nrow(check_data) > 1){
        input_data$state_from[which(input_data$prob_matrix == "initial")] <- apply(initialProbs, 1, draw_sample, missing_variable_levels = missing_variable_levels, without_trans_prob = without_trans_prob)
        input_data$prob_matrix[which(input_data$prob_matrix == "initial")] <- "forward"

        for (iYear in 1:number_of_transitions){


          if (is.null(covariates_transition)){
            newdata = data.frame(state_from = input_data[which(input_data$prob_matrix == "forward"& input_data$tran_Year == iYear & is.na(input_data$state_from) == F), c("state_from")])
          } else {
            newdata = input_data[which(input_data$prob_matrix == "forward"& input_data$tran_Year == iYear & is.na(input_data$state_from) == F), c("state_from", covariates_transition)]
          }

          forwardProbabilities <- predict(modelFittransition, newdata = newdata, "prob")
          input_data$state_to[which(input_data$prob_matrix == "forward" & input_data$tran_Year == iYear & is.na(input_data$state_from) == F)] <- apply(forwardProbabilities, 1, draw_sample, missing_variable_levels = missing_variable_levels, without_trans_prob = without_trans_prob)
          if (iYear == nrow(input_data)) break
          input_data$state_from[which(input_data$tran_Year == iYear+1)] <- input_data$state_to[which(input_data$tran_Year == iYear)]

        }
      }

      if (nrow(check_data) == 1){

        input_data$state_from[which(input_data$prob_matrix == "initial")] <- draw_sample(initialProbs, missing_variable_levels = missing_variable_levels, without_trans_prob = without_trans_prob)
        input_data$prob_matrix[which(input_data$prob_matrix == "initial")] <- "forward"

        for (iYear in 1:number_of_transitions){
          if (iYear == 1){

            if (is.null(covariates_transition)){
              newdata = data.frame(state_from = input_data[which(input_data$prob_matrix == "forward"& input_data$tran_Year == iYear & is.na(input_data$state_from) == F), c("state_from")])
            } else {
              newdata = input_data[which(input_data$prob_matrix == "forward"& input_data$tran_Year == iYear & is.na(input_data$state_from) == F), c("state_from", covariates_transition)]
            }

            if (nrow(newdata) == 1){

              forwardProbabilities <- predict(modelFittransition, newdata = newdata, "prob")
              input_data$state_to[which(input_data$prob_matrix == "forward" & input_data$tran_Year == iYear & is.na(input_data$state_from) == F)] <- draw_sample(forwardProbabilities, missing_variable_levels = missing_variable_levels, without_trans_prob = without_trans_prob)
              if (iYear == nrow(input_data)) break
              input_data$state_from[which(input_data$tran_Year == iYear+1)] <- input_data$state_to[which(input_data$tran_Year == iYear)]
            }
            if (nrow(newdata) > 1){

              forwardProbabilities <- predict(modelFittransition, newdata = newdata, "prob")
              input_data$state_to[which(input_data$prob_matrix == "forward" & input_data$tran_Year == iYear & is.na(input_data$state_from) == F)] <- apply(forwardProbabilities, 1, draw_sample, missing_variable_levels = missing_variable_levels, without_trans_prob = without_trans_prob) #draw_sample(forwardProbabilities, missing_variable_levels)
              if (iYear == nrow(input_data)) break
              input_data$state_from[which(input_data$tran_Year == iYear+1)] <- input_data$state_to[which(input_data$tran_Year == iYear)]
            }
          }
          if (iYear > 1){

            if (is.null(covariates_transition)){
              newdata = data.frame(state_from = input_data[which(input_data$prob_matrix == "forward"& input_data$tran_Year == iYear & is.na(input_data$state_from) == F), c("state_from")])
            } else {
              newdata = input_data[which(input_data$prob_matrix == "forward"& input_data$tran_Year == iYear & is.na(input_data$state_from) == F), c("state_from", covariates_transition)]
            }

            if (nrow(newdata) == 1){

              forwardProbabilities <- predict(modelFittransition, newdata = newdata, "prob")
              input_data$state_to[which(input_data$prob_matrix == "forward" & input_data$tran_Year == iYear & is.na(input_data$state_from) == F)] <- draw_sample(forwardProbabilities, missing_variable_levels)
              if (iYear == nrow(input_data)) break
              input_data$state_from[which(input_data$tran_Year == iYear+1)] <- input_data$state_to[which(input_data$tran_Year == iYear)]
            }
            if (nrow(newdata) > 1){

              forwardProbabilities <- predict(modelFittransition, newdata = newdata, "prob")
              input_data$state_to[which(input_data$prob_matrix == "forward" & input_data$tran_Year == iYear & is.na(input_data$state_from) == F)] <- apply(forwardProbabilities, 1, draw_sample, missing_variable_levels = missing_variable_levels, without_trans_prob = without_trans_prob) #draw_sample(forwardProbabilities, missing_variable_levels)
              if (iYear == nrow(input_data)) break
              input_data$state_from[which(input_data$tran_Year == iYear+1)] <- input_data$state_to[which(input_data$tran_Year == iYear)]
            }

          }
        }
      }
      return(input_data)

    }

    check_data <- input_data[which(input_data$prob_matrix == "initial"), covariates_initial]
    if (nrow(check_data) == 0){
      initial_forward_data <- input_data
    } else {
      initial_forward_data <- initial_forward_function(input_data = input_data,
                                                       number_of_transitions = number_of_transitions,
                                                       covariates_initial = covariates_initial,
                                                       covariates_transition = covariates_transition,
                                                       missing_variable_levels = missing_variable_levels)
    }
  } else {


    initial_forward_function <- function(input_data = input_data,
                                         number_of_transitions = number_of_transitions,
                                         covariates_initial = covariates_initial,
                                         covariates_transition = covariates_transition,
                                         missing_variable_levels = missing_variable_levels){


      for (iYear in 1:number_of_transitions){
        if (iYear == 1){

          if (is.null(covariates_transition)){
            newdata = data.frame(state_from = input_data[which(input_data$prob_matrix == "forward"& input_data$tran_Year == iYear & is.na(input_data$state_from) == F), c("state_from")])
          } else {
            newdata = input_data[which(input_data$prob_matrix == "forward"& input_data$tran_Year == iYear & is.na(input_data$state_from) == F), c("state_from", covariates_transition)]
          }
          if (nrow(newdata) == 1){

            forwardProbabilities <- predict(modelFittransition, newdata = newdata, "prob")
            input_data$state_to[which(input_data$prob_matrix == "forward" & input_data$tran_Year == iYear & is.na(input_data$state_from) == F)] <- draw_sample(forwardProbabilities, missing_variable_levels)
            if (iYear == nrow(input_data)) break
            input_data$state_from[which(input_data$tran_Year == iYear+1)] <- input_data$state_to[which(input_data$tran_Year == iYear)]
          }
          if (nrow(newdata) > 1){

            forwardProbabilities <- predict(modelFittransition, newdata = newdata, "prob")
            input_data$state_to[which(input_data$prob_matrix == "forward" & input_data$tran_Year == iYear & is.na(input_data$state_from) == F)] <- apply(forwardProbabilities, 1, draw_sample, missing_variable_levels = missing_variable_levels) #draw_sample(forwardProbabilities, missing_variable_levels)
            if (iYear == nrow(input_data)) break
            input_data$state_from[which(input_data$tran_Year == iYear+1)] <- input_data$state_to[which(input_data$tran_Year == iYear)]
          }
        }
        if (iYear > 1){

          if (is.null(covariates_transition)){
            newdata = data.frame(state_from = input_data[which(input_data$prob_matrix == "forward"& input_data$tran_Year == iYear & is.na(input_data$state_from) == F), c("state_from")])
          } else {
            newdata = input_data[which(input_data$prob_matrix == "forward"& input_data$tran_Year == iYear & is.na(input_data$state_from) == F), c("state_from", covariates_transition)]
          }

          if (nrow(newdata) == 1){

            forwardProbabilities <- predict(modelFittransition, newdata = newdata, "prob")
            input_data$state_to[which(input_data$prob_matrix == "forward" & input_data$tran_Year == iYear & is.na(input_data$state_from) == F)] <- draw_sample(forwardProbabilities, missing_variable_levels = missing_variable_levels, without_trans_prob = without_trans_prob)
            if (iYear == nrow(input_data)) break
            input_data$state_from[which(input_data$tran_Year == iYear+1)] <- input_data$state_to[which(input_data$tran_Year == iYear)]
          }
          if (nrow(newdata) > 1){

            forwardProbabilities <- predict(modelFittransition, newdata = newdata, "prob")
            input_data$state_to[which(input_data$prob_matrix == "forward" & input_data$tran_Year == iYear & is.na(input_data$state_from) == F)] <- apply(forwardProbabilities, 1, draw_sample, missing_variable_levels = missing_variable_levels, without_trans_prob = without_trans_prob) #draw_sample(forwardProbabilities, missing_variable_levels)
            if (iYear == nrow(input_data)) break
            input_data$state_from[which(input_data$tran_Year == iYear+1)] <- input_data$state_to[which(input_data$tran_Year == iYear)]
          }
        }
      }
      return(input_data)

    }
  }
  initial_forward_data <- initial_forward_function(input_data = input_data,
                                                   number_of_transitions = number_of_transitions,
                                                   covariates_initial = covariates_initial,
                                                   covariates_transition = covariates_transition,
                                                   missing_variable_levels = missing_variable_levels)




  #### 2.2 intermittent function ####
  imputeIntermittent <- function(input_data = input_data,
                                 patient_id = patient_id,
                                 from_state = from_state,
                                 to_state = to_state,
                                 year_missing = year_missing,
                                 covariates_initial = covariates_initial,
                                 covariates_transition = covariates_transition,
                                 missing_variable_levels = missing_variable_levels){


    # (A1) Create the function extratProb_new, which is amazing (change name "extratProb_new", and also explain why it is amazing)
    extratProb_new <- function(diff_scenarios, prob_list){
      trans_probs <- c()
      trans_probs[1] <- prob_list[[1]][diff_scenarios[1,2]]
      trans_probs <- unlist(trans_probs)
      for (iROW in 2:length(prob_list)){
        trans_probs[iROW] <- prob_list[[iROW]][diff_scenarios[iROW,1],diff_scenarios[iROW,2]]
      }
      return(prod(trans_probs))
    }

    # (A2) Need function to calculate all scenarios #
    computeIndexTrans <- function(from_state, to_state, year_missing){
      allPossibilities <- expand.grid(rep(list(1:length(missing_variable_levels)),year_missing))
      finalScenarios <- as.matrix(cbind(rep(from_state, nrow(allPossibilities)), allPossibilities, rep(to_state, nrow(allPossibilities))))
      colnames(finalScenarios) <- rownames(finalScenarios) <- NULL
      return(finalScenarios)
    }

    # (A3) Need function to change Word to Numeric function #
    changeWord2Num <- function(word2change, missing_variable_levels){
      return(which(missing_variable_levels == word2change))
    }

    ## computing the probability of each senario
    from_state <- changeWord2Num(from_state, missing_variable_levels)
    to_state <- changeWord2Num(to_state, missing_variable_levels)
    allScenarios <- computeIndexTrans(from_state, to_state, year_missing)
    probs_all_scenarios <- rep(NA, nrow(allScenarios))



    prob_list <- list()
    if (is.null(covariates_transition)){
      newdata = data.frame(state_from = input_data[1, c("state_from")])
    } else {
      newdata = input_data[1, c("state_from", covariates_transition)]
    }
    prob_list[[1]] <- predict(modelFittransition, newdata = newdata, "prob")

    for (iRow in 2:nrow(input_data)){

      data2pred_list <- list()
      for (ilevel in 1:length(missing_variable_levels)){
        data2pred_list[[ilevel]] <- data.frame(state_from = missing_variable_levels[ilevel], input_data[iRow, covariates_transition])
      }
      data2pred <- do.call("rbind", data2pred_list)

      prob_list[[iRow]] <- predict(modelFittransition, newdata = data2pred, "prob")

    }
    prob_list


    probs_all_scenarios <- c()

    for (iRow in 1:nrow(allScenarios)){
      tmpScen <- allScenarios[iRow,]
      allidx <- cbind(tmpScen[1:(length(tmpScen)-1)], tmpScen[2:length(tmpScen)])
      probs_all_scenarios[iRow] <- extratProb_new(allidx, prob_list)
    }

    normalizedProbs <- probs_all_scenarios/sum(probs_all_scenarios)
    imputedScenario <- sample(1:nrow(allScenarios), 1, prob = normalizedProbs)
    imp2return <- allScenarios[imputedScenario,]
    imp2returnInWords <- rep(NA,length(imp2return))

    find_levels = unique(imp2return)
    for  (iWords in 1:length(find_levels)){
      imp2returnInWords[which(imp2return == find_levels[iWords])] <- missing_variable_levels[find_levels[iWords]]
    }
    ## These lines are for making scenarios in words if user is interested
    # allScenariosInWords <- allScenarios
    # allScenariosInWords[which(allScenarios == 1)] <- missing_variable_levels[1]
    # allScenariosInWords[which(allScenarios == 2)] <- missing_variable_levels[2]
    # allScenariosInWords[which(allScenarios == 3)] <- missing_variable_levels[3]
    # allScenariosInWords[which(allScenarios == 4)] <- missing_variable_levels[4]
    # allScenariosInWords[which(allScenarios == 5)] <- missing_variable_levels[5]
    # allScenariosInWords[which(allScenarios == 6)] <- missing_variable_levels[6]
    # allScenariosInWords[which(allScenarios == 7)] <- missing_variable_levels[7]
    # allScenariosInWords[which(allScenarios == 8)] <- missing_variable_levels[8]
    # allScenariosInWords[which(allScenarios == 9)] <- missing_variable_levels[9]
    # allScenariosInWords[which(allScenarios == 10)] <- missing_variable_levels[10]
    # rownames(allScenariosInWords) <- paste("Scenario", 1:nrow(allScenariosInWords))
    # names(normalizedProbs) <- paste("Scenario", 1:nrow(allScenariosInWords))
    #return(list(scenarios  = allScenariosInWords, probs  = normalizedProbs, imputedData = imp2returnInWords[-c(1,length(imp2returnInWords))]))

    imputedData = imp2returnInWords[-c(1,length(imp2returnInWords))]
    for (iRow in 1:nrow(input_data)){

      if (is.na(input_data$state_to[iRow] == T)){
        input_data$state_to[iRow] <- imputedData[iRow]
      }
      if (iRow == nrow(input_data)) break
      input_data$state_from[iRow+1] <- imputedData[iRow]
    }
    return(input_data)
  }

  inter_missing <- initial_forward_data
  inter_missing <- inter_missing[inter_missing$prob_matrix == "intermittent",]

  if (nrow(inter_missing) > 0){
    inter_missing$id_inter <- rep(0, nrow(inter_missing))
    inter_missing$id_inter[which(!is.na(inter_missing$state_from))] <- 1
    inter_missing$id_inter = cumsum(inter_missing$id_inter)

    idx.subj = unique(inter_missing[, patient_id])
    bdd.changed_inter = inter_missing
    for (i in 1:length(idx.subj)){
      input.subj <- bdd.changed_inter[bdd.changed_inter[, patient_id] == idx.subj[i],]
      inter_data = imputeIntermittent(input_data = input.subj, patient_id = patient_id, input.subj$state_from[1], input.subj$state_to[nrow(input.subj)], nrow(input.subj)-1,
                                      covariates_initial = covariates_initial,
                                      covariates_transition = covariates_transition,
                                      missing_variable_levels = missing_variable_levels)
      bdd.changed_inter[bdd.changed_inter[, patient_id] == idx.subj[i],] = inter_data
    }
    bdd.changed_inter$id_inter <- NULL

    initial_forward_data <- initial_forward_data[initial_forward_data$prob_matrix != "intermittent",]
    input_data <- rbind(initial_forward_data, bdd.changed_inter)

  } else {
    input_data <- initial_forward_data
  }


  #### 2.3 backward function ####
  backward_function <- function(input_data = input_data,
                                patient_id = patient_id,
                                number_of_transitions = number_of_transitions,
                                covariates_initial = covariates_initial,
                                covariates_transition = covariates_transition,
                                missing_variable_levels = missing_variable_levels){

    for (iYear in number_of_transitions:1){

      # Get initialProbs for each tran_Year starting from the last (i.e. 1st time that backward observation occurred)
      initialProbs <- predict(modelFitinitial, newdata = input_data[which(input_data$tran_Year == iYear), covariates_initial], "prob")

      # Create pseudo-observations for this patient for each tran_Year
      data2pred_list <- list()
      for (ilevel in 1:length(missing_variable_levels)){
        data2pred_list[[ilevel]] <- data.frame(state_from = missing_variable_levels[ilevel], input_data[which(input_data$tran_Year == iYear), covariates_transition])
      }
      data2pred <- do.call("rbind", data2pred_list)

      # Calculate transitionProbs on those pseudo-data
      transitionProbs <- predict(modelFittransition, newdata = data2pred, "prob")
      rownames(transitionProbs) <- missing_variable_levels

      # Compute joint probabilities as product of (transitionProbs, initialProbs)
      jointProbabilities_list <- list()
      for (ilevel in 1:length(missing_variable_levels)){
        jointProbabilities_list[[ilevel]] <-  transitionProbs[ilevel,] * initialProbs[[ilevel]]
      }
      jointProbabilities <-  do.call("rbind", jointProbabilities_list)

      # Compute state_to initial probabilities
      initialProbsStateTo <- apply(jointProbabilities, 2, sum)

      #  Compute backward probabilities as product of (transitionProbs, initialProbs)
      backwardProbabilities_list <- list()
      for (ilevel in 1:length(missing_variable_levels)){
        backwardProbabilities_list[[ilevel]] <-  jointProbabilities[,ilevel] / initialProbsStateTo[ilevel]
      }
      backwardProbabilities <-  do.call("cbind", backwardProbabilities_list)

      colnames(backwardProbabilities) <- missing_variable_levels; backwardProbabilities

      status <- input_data$state_to[which(input_data$tran_Year == iYear)]

      probs_back = backwardProbabilities[, status]

      input_data$state_from[which(input_data$tran_Year == iYear)] <- draw_sample(probs_back, missing_variable_levels = missing_variable_levels, without_trans_prob = without_trans_prob)
      input_data$state_to[which(input_data$tran_Year == iYear-1)] <- input_data$state_from[which(input_data$tran_Year == iYear)]

      # } else  {
      #   warning(paste("For patient ", patient_id, " equal probabilities were used to generate backward probabilities"), sep = "")
      #   input_data$state_from[which(input_data$tran_Year == iYear)] <- sample(missing_variable_levels, 1)
      #   input_data$state_to[which(input_data$tran_Year == iYear-1)] <- input_data$state_from[which(input_data$tran_Year == iYear)]
      #
      # }


    }
    return(input_data)
  }

  back_missing <- initial_forward_data
  back_missing <- back_missing[back_missing$prob_matrix == "backward",]

  if (nrow(back_missing) > 0){
    idx.subj = unique(back_missing[, patient_id])

    bdd.changed_back = back_missing

    for (i in 1:length(idx.subj)){
      input.subj <- bdd.changed_back[bdd.changed_back[, patient_id] == idx.subj[i],]
      number_of_transitions <- nrow(input.subj)
      back_data = backward_function(input_data = input.subj, patient_id = patient_id, number_of_transitions = number_of_transitions,
                                    covariates_initial = covariates_initial,
                                    covariates_transition = covariates_transition,
                                    missing_variable_levels = missing_variable_levels)
      bdd.changed_back[bdd.changed_back[, patient_id] == idx.subj[i],] = back_data
    }


    input_data <- input_data[input_data$prob_matrix != "backward",]
    input_data <- rbind(input_data, bdd.changed_back)
    input_data <- input_data[order(input_data[, patient_id], input_data$tran_Year),]
  } else {
    input_data <- input_data
  }
  return(list(input_data= input_data, modelFittransition= modelFittransition, modelFitinitial = modelFitinitial,
              observed_smoke_transition = observed_smoke_transition,
              observed_smoke_initial = observed_smoke_initial))

}


















