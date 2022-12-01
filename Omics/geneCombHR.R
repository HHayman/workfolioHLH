library(survival)


geneCombHR <-
  function(maf = NULL,
           genes = NULL,
           maxN = 5,
           Time = NULL,
           Event = NULL,
           CImin = 0.75,
           CImax = 2,
           minMut = 5,
           adjust = "No",
           fileName = NULL,
           filePath = NULL) {
    ## Extract clinical and maf data

    mafData <- maf@data
    clinData <- maf@clinical.data

    if (is.null(genes)) {
      genes = getGeneSummary(x = maf)[1:20, Hugo_Symbol]
    }

    clinData$Time = suppressWarnings(as.numeric(as.character(clinData[[Time]])))
    clinData$Event = suppressWarnings(as.integer(clinData[[Event]]))
    clinData$Time = ifelse(test = is.infinite(clinData$Time),
                           yes = NA,
                           no = clinData$Time)

    genes <<- genes

    ## Define gene combinations
    testResults <- data.frame()

    for (i in 1:maxN) {
      geneCombs <- data.frame(combn(x = genes, m = i))

      ## Run coxph
      groupNames = c("Mutant", "WT")

      zphP <- vector()
      zphComb <- vector()

      for (j in 1:ncol(geneCombs)) {
        coxGenes <- geneCombs[, j]

        vecTSB <- vector()
        for (k in 1:length(coxGenes)) {
          tempTSB <-
            as.character(unique(mafData[mafData$Hugo_Symbol %in% coxGenes[k], ]$Tumor_Sample_Barcode))
          vecTSB <- c(vecTSB, tempTSB)
        }
        v1 <- table(vecTSB)
        coxTSB <- names(v1)[v1 == k]

        clinData$Group <-
          ifelse(clinData$Tumor_Sample_Barcode %in% coxTSB,
                 groupNames[1],
                 groupNames[2])

        clinData$Group <-
          factor(clinData$Group, levels = c("WT", "Mutant"))

        clinData <<- clinData

        if (nrow(clinData[clinData$Group == "Mutant", ]) >= minMut) {
          surv.cox <-
            survival::coxph(formula = survival::Surv(time = Time, event = Event) ~ Group,
                            data = clinData)

          surv.cox.test <- cox.zph(surv.cox)
          zphP[[j]] <- surv.cox.test$table[ , "p" ][[2]]
          zphComb[[j]] <- paste(coxGenes, collapse = "_")

          tempVec <- vector()
          tempVec[1] <- paste(coxGenes, collapse = "_")
          tempVec[2] <-
            round(summary(surv.cox)$coefficients[, 5], 3)
          tempVec[3] <- round(exp(coef(surv.cox)), 2)
          tempVec[4] <- round(exp(confint(surv.cox))[1], 2)
          tempVec[5] <- round(exp(confint(surv.cox))[2], 2)
          tempVec[6] <- nrow(clinData[clinData$Group == "WT"])
          tempVec[7] <- nrow(clinData[clinData$Group == "Mutant"])

          testResults <- rbind(testResults, t(tempVec))
        }
      }
    }

    testResults <<- testResults

    colnames(testResults) <-
      c("Genes_combinations",
        "P_value",
        "HR",
        "lowCI",
        "highCI",
        "WT",
        "Mutant")
    testResults$adj.p <- p.adjust(testResults$P_value, method = "fdr")
    testResults <<- testResults

    Filepath <- paste(filePath, fileName, "testResults.csv", sep = "-")
    write.csv(testResults, Filepath)


    sigTestResults <- testResults
    sigTestResults$hrCode <-
      ifelse(sigTestResults$HR < 0.5 | sigTestResults$HR > 2, 1, 0)
    sigTestResults$pCode <-
      ifelse(sigTestResults$P_value < 0.05, 1, 0)
    sigTestResults <-
      sigTestResults[sigTestResults$hrCode == 1 &
                       sigTestResults$pCode == 1,]
    sigTestResults <<- sigTestResults


    sigAdjResults <- testResults[testResults$adj.p < 0.05,]

    Filepath <- paste(filePath, fileName, "sigTestResults.csv", sep = "-")
    write.csv(sigTestResults, Filepath)

    zphData <- data.frame(zphComb, zphP)

    setClass("combHR",
             representation(Results = "data.frame", sigResults = "data.frame", sigAdjResults = "data.frame", cox.zph.res = "data.frame"))
    output <-
      new("combHR", Results = testResults, sigResults = sigTestResults, sigAdjResults = sigAdjResults, cox.zph.res = zphData)

    return(output)

  }



