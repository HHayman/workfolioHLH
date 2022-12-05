

library(dplyr)
library(magrittr)
library(maftools)




geneStatus <-
  function(maf = NULL,
           ID = NULL,
           genes = NULL) {
    
    variantClassificationGenesTotal <- maf@gene.summary %>% arrange(Hugo_Symbol)
    
    #### Variant classifications parent data
    variantsTable <- maf@data %>% select(any_of(c(ID, "Hugo_Symbol", "Variant_Classification")))
    dataFramesParent <-
      setNames(
        data.frame(matrix(ncol = 8, nrow = 0)),
        c(
          "Gene",
          "Frame_Shift_Del",
          "Frame_Shift_Ins",
          "In_Frame_Del",
          "In_Frame_Ins",
          "Missense_Mutation",
          "Nonsense_Mutation",
          "Splice_Site"
        )
      )
    IDs <- unique(variantsTable[[ID]]) %>% lapply(as.character)
    for (i in 1:length(IDs)) {
      caseData <- variantsTable %>% filter(!!sym(ID) == IDs[i])
      Genes <- unique(caseData$Hugo_Symbol)
      for (j in 1:length(Genes)) {
        geneData <- caseData %>% filter(Hugo_Symbol == Genes[j]) %>% select(Variant_Classification) %>% table() %>% data.frame()
        colnames <- geneData$Variant_Classification
        geneData <- t(geneData[,-1]) %>% `colnames<-`(colnames)
        geneData <-
          cbind(Gene = Genes[j], ID = IDs[i], geneData)
        dataFramesParent <- rbind(dataFramesParent, geneData)
      }
    }
    
    
    # Get row number
    rows <- nrow(dataFramesParent)
    
    # Create vector to hold codes
    zeroCount <- vector()
    
    # Generate multihit codes
    for (i in 1:nrow(dataFramesParent)) {
      tempData <- dataFramesParent[i,] %>% t() %>% data.frame()
      tabRow <- tempData[3:9,] %>% table() %>% data.frame()
      ifelse(tabRow[1,]$Freq < 6, zeroCount[i] <-
               1, zeroCount[i] <- 0)
    }
    dataFramesParent$multiCode <- zeroCount
    
    
    # Set individual codes
    dataFramesParent <- dataFramesParent %>%
      mutate(indivCode = 0) %>%
      mutate(
        indivCode = case_when(
          Frame_Shift_Del > 0 ~ 1,
          Frame_Shift_Ins > 0 ~ 2,
          In_Frame_Del > 0 ~ 3,
          In_Frame_Ins > 0 ~ 4,
          Missense_Mutation > 0 ~ 5,
          Nonsense_Mutation > 0 ~ 6,
          Splice_Site > 0 ~ 7
        )
      )
    
    
    # Set parent code
    dataFramesParent$variantCode <-
      ifelse(dataFramesParent$multiCode == 1,
             8,
             dataFramesParent$indivCode)
    
    # Arrange columns
    dataFramesParent <- dataFramesParent %>% select(-c(multiCode, indivCode)) %>% dplyr::relocate(variantCode, .after = ID)
    
    ## Summing scores
    dataFramesParent[, 4:10] <-
      lapply(dataFramesParent[, 4:10], as.numeric)
    dataFramesParent <-
      dataFramesParent %>%
      mutate(TotalMut = rowSums(dataFramesParent[, 4:10])) %>%
      filter(Gene %in% genes) %>%
      mutate(Gene = as.character(Gene), ID = as.character(ID)) %>%
      data.frame()
    
    
    
    return(dataFramesParent)
    
    
    
    
    
  }
