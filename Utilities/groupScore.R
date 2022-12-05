library(magrittr)
library(dplyr)



groupScore <-
  function(data = NULL,
           groupNames = NULL,
           groupVars = NULL,
           groupValidator = NULL,
           groupValidatorValues = NULL,
           scoreName = NULL
           ) {

    breakValue <- round(1 / length(unlist(groupVars)), 2)

    varsA <- unlist(groupVars[[1]])
    varsB <- unlist(groupVars[[2]])
    vars <- c(varsA, varsB)

    data <- data %>%
      mutate_at(varsA, funs(ifelse(. == 1, breakValue, .))) %>%
      mutate_at(varsB, funs(ifelse(. == 1, breakValue * -1, .))) %>%
      dplyr::mutate(score = rowSums(across(vars)), scorePred = ifelse(score > 0, groupNames[1], ifelse(score < 0, groupNames[2], NA_real_))) %>%
      rename_at(vars(score), funs(paste0(scoreName))) %>%
      rename_at(vars(scorePred), funs(paste0(scoreName, "Pred", sep = "")))

    return(data)

  }
