#' #'Survflow
#' 
#' 
#' #' @author Hannah Hayman, \email{hannah.louise.hayman@@gmail.com}
#' #' This function takes in a dataset(s) and runs appropriate survival analysis.
#' 
#' 
#' #' @import survival
#' #' @import maxstat
#' #' @import survminer
#' #' @import tidyverse
#' #' @import ggpubr
#' #' @import ggplot2
#' #' @import gplots
#' #' @import flextable
#' 
#' 
#' #' @param Data A Coded dataset
#' #' @param Variables A list of variables (coded as 0 and 1) for analysis
#' #' @param LegendLabels Optional labels for legend
#' #' @param Identifier Identifier variable for cases
#' #' @param PlotTitles Optional plot titles
#' #' @param SurvivalStatus A status variable (coded as 0 and 1) for survival analysis
#' #' @param SurvivalTime A survival time variable - continuous
#' #' @param SurvivalTimeUnit The unit of time for survival time. Months etc
#' #' @param xYearSurvivalVar The number of years to be used to calculate 'X' years survival. Default = 5
#' #' @param SurvBase Takes TRUE or FALSE. Toggle to activate base survival analysis.
#' 
#' 
#' #' @return Takes a coded dataset (use CutFlow) and runs survival analysis for each variable, by each survival outcome, and returns output.
#' 
#' 
#' #' @examples SurvFlow(Glasgow,  Variables = c("GD_PercPositiveCellsinHealthyTissue_Coded", "GD_PercPositiveCellsinTumourTissue_Coded", "CD8_PercPositiveCellsinHealthyTissue_Coded", "CD8_PercPositiveCellsinTumourTissue_Coded"), LegendLabels = c("Low", "High"), Identifier = "TMA_ID", PlotTitles = c("GD Healthy Stroma", "GD Healthy Epithelium", "GD Healthy Tissue", "GD Tumour Stroma", "GD Tumour Epithelium", "GD Tumour Tissue", "CD8 Healthy Stroma", "CD8 Healthy Epithelium", "CD8 Healthy Tissue", "CD8 Tumour Stroma", "CD8 Tumour Epithelium", "CD8 Tumour Tissue"), SurvivalStatus = c("CSS", "OS", "DFS", "RFS"), SurvivalTime = c("CSS_2017", "CSS_2017", "DFSmonths", "DFSmonths"), xYearSurvivalVar = 5, SurvivalTimeUnit = "Months", SurvBase = TRUE)
#' 
#' #' @export



# Packages
library(survival)
library(maxstat)
library(survminer)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(gplots)
library(grid)
library(gridExtra)
library(flextable)


# Function
SurvFlow <-
  function(Data,
           Variables,
           LegendLabels = NULL,
           Identifier,
           PlotTitles = NULL,
           SurvivalStatus,
           SurvivalTime,
           SurvivalTimeUnit,
           xYearSurvivalVar = 5,
           SurvBase = FALSE) {
    #### Read data and check for errors ####

    # Save original working directory and set to restore on exit
    WD <- toString(getwd())
    on.exit(setwd(WD), add = TRUE)

    # Check that dataset is a dataframe
    if (!inherits(Data, "data.frame"))
      #stop("data should be an object of class data.frame")
      Data <- as.data.frame(Data)

    # Check that identifier is present in the dataset
    if (!all(Identifier %in% colnames(Data)))
      stop("Identifier not found in the data: ")

    # Check that survival variables are in the dataset
    if (!all(c(SurvivalStatus, SurvivalTime) %in% colnames(Data)))
      stop("Specify correct column names containing SurvivalStatus and SurvivalTime values.")

    # Check that variable names are present in the dataset
    if (!all(Variables %in% colnames(Data)))
      stop("Some variables are not found in the data: ",
           paste(setdiff(Variables, colnames(Data)), collapse = ", "))



    #### Create a subdirectory for output ####

    # Define and check subdirectory name
    DataFileName <- deparse(substitute(Data))
    Number <- 1
    OutputDirectory <-
      paste("SurvFlow", DataFileName, Sys.Date(), Number, sep = '_')
    while (file.exists(OutputDirectory))
    {
      Number <- Number + 1
      OutputDirectory <-
        paste("SurvFlow", DataFileName, Sys.Date(), Number, sep = '_')
    }

    # Create the subdirectory
    dir.create(OutputDirectory)
    setwd(OutputDirectory)
    OD <- toString(getwd())




    #### Baseline survival analysis, per status, per variable ####

    # Check if base survival is requested
    if (SurvBase == TRUE) {
      # Create and set a subdirectory
      BaseSurvDirectoryName <-
        paste("BaseSurvival", sep = "")
      dir.create(BaseSurvDirectoryName)
      setwd(BaseSurvDirectoryName)
      BSD <- toString(getwd())

      ## Median table
      medianData <- data.frame()
      zphDataList <- list()
      zphP <- vector()
      zphSurv <- vector()
      zphVar <- vector()

      for (j in 1:length(SurvivalStatus)) {
        # Create and set a subdirectory
        DirectoryName <-
          paste("BaseSurvival_", SurvivalStatus[j], sep = "")
        dir.create(DirectoryName)
        setwd(DirectoryName)


        for (i in 1:length(Variables)) {
          #  Create temporary dataset
          TempData <-
            subset(Data,
                   select = c(
                     Identifier,
                     SurvivalStatus[j],
                     SurvivalTime[j],
                     Variables[i]
                   ))

          TimeName <- SurvivalTime[j]
          StatusName <- SurvivalStatus[j]
          IDName <- Identifier
          VarName <- Variables[i]

          colnames(TempData)[4] <- "Variable"
          colnames(TempData)[2] <- "Status"
          colnames(TempData)[3] <- "Time"


          TempData <- na.omit(TempData)

          groupCheck <- as.data.frame(table(unique(TempData$Variable)))


          # if (nrow(groupCheck) == 2){
          # Events <- TempData %>% group_by(Variable) %>%
          #   summarise(Events_Sum = sum(Status))
          #
          # Cases <- TempData %>% group_by(Variable) %>% tally()



          if (nrow(groupCheck) == 2){
            Events <- TempData %>% dplyr::group_by(Variable) %>%
              dplyr::summarise(Events_Sum = sum(Status))

            Cases <- TempData %>% dplyr::group_by(Variable) %>% tally()


          SurvOutputTable <-
            merge(Cases, Events, by.x = "Variable", all.y = TRUE)
          colnames(SurvOutputTable)[1] <- "Group"
          colnames(SurvOutputTable)[2] <- "Cases"
          colnames(SurvOutputTable)[3] <- "Events"

          LabX <-
            paste(StatusName, " (", SurvivalTimeUnit, ")", sep = "")
          LabY <- "Survival Outcome Probability"
          ifelse (
            is.null(LegendLabels),
            LegLabs <- unique(TempData$Variable),
            LegLabs <- LegendLabels
          )

          if (nrow(SurvOutputTable) == 2) {
            ## Calculate X-year survival and add to SurvOutputTable table
            # Create independent tables for each group
            TempData0 <- subset(TempData, Variable == 0)
            TempData1 <- subset(TempData, Variable == 1)

            # Calculate times value
            xTime <- 12 * xYearSurvivalVar

            # Fit survival by x # of years
              xYearSurvival <-
                summary(survfit(Surv(Time, Status) ~ Variable, data = TempData),
                        times = xTime,
                        extend = TRUE)
              xYearSurvival0 <-
                summary(survfit(Surv(Time, Status) ~ Variable, data = TempData0),
                        times = xTime,
                        extend = TRUE)
              xYearSurvival1 <-
                summary(survfit(Surv(Time, Status) ~ Variable, data = TempData1),
                        times = xTime,
                        extend = TRUE)

            # Define x-year survival values
            xYearSurvivalValues <-
              c(
                paste(
                  (round(xYearSurvival0$surv, digits = 2) * 100),
                  "% (",
                  (round(xYearSurvival0$lower, digits = 2) * 100),
                  "%, ",
                  (round(xYearSurvival0$upper, digits = 2) * 100),
                  "%)",
                  sep = ""
                ),
                paste(
                  (round(xYearSurvival1$surv, digits = 2) * 100),
                  "% (",
                  (round(xYearSurvival1$lower, digits = 2) * 100),
                  "%, ",
                  (round(xYearSurvival1$upper, digits = 2) * 100),
                  "%)",
                  sep = ""
                )
              )

            # Define variable name for x-year survival and add to SurvOutputTable table
            xYearSurvivalName <-
              paste(xYearSurvivalVar, "-Year Survival", sep = "")
            SurvOutputTable$xYearSurvivalColumn <-
              xYearSurvivalValues
            colnames(SurvOutputTable)[4] <- xYearSurvivalName
          }


          ## Run cox hazard analysis and add to SurvOutputTable table
          # Fit univariate cox hazard regression
          UniCoxFit <-
            coxph(Surv(Time, Status) ~ Variable, data = TempData) %>% gtsummary::tbl_regression(exp = TRUE)
          UniCoxFitTable <- UniCoxFit$table_body


          ## zph Tests
          surv.cox.test <- cox.zph(coxph(Surv(Time, Status) ~ Variable, data = TempData))
          zphP[[i]] <- surv.cox.test$table[ , "p" ][[2]]
          zphSurv[[i]] <- SurvivalStatus[j]
          zphVar[[i]] <- Variables[i]




          # Add univariate cox regression to table
          UniCoxValuesHR <-
            paste(
              round(UniCoxFitTable$estimate, digits = 2),
              " (",
              round(UniCoxFitTable$conf.low, digits = 2),
              ", ",
              round(UniCoxFitTable$conf.high, digits = 2),
              ")",
              sep = ""
            )
          UniCoxValuesP <-
            paste(round(UniCoxFitTable$p.value, digits = 2), sep = "")

          SurvOutputTable$Cox_HR <- UniCoxValuesHR
          SurvOutputTable$Cox_HR_P <- UniCoxValuesP



          # Fit multivariate cox hazard regression



          # Calculate logrank and add to SurvOutputTable table
          LogRankData <-
            survdiff(Surv(Time, Status) ~ Variable, data = TempData)
          LogRankP <-
            pchisq(LogRankData$chisq,
                   length(LogRankData$n) - 1,
                   lower.tail = FALSE)
          SurvOutputTable$Logrank_P <- round(LogRankP, digits = 3)


          ## Additional formatting of SurvOutputTable table
          # Rename the groups in SurvOutputTable table
          SurvOutputTable$Group <- LegLabs



          ### Create median survival table
          tempVec <- vector()

          fit <- survfit(Surv(Time, Status) ~ Variable, data = TempData)
          medSurv <- surv_median(fit)



          tempVec[1] <- substr(deparse(substitute(Data)), 1, nchar(deparse(substitute(Data))) -5)
          tempVec[2] <- SurvivalStatus[j]
          tempVec[3] <- substr(Variables[i], 1, nchar(Variables[i]) -6)
          tempVec[4] <- format(round(medSurv[1, 2], 2), nsmall = 2)
          tempVec[5] <- format(round(medSurv[2, 2], 2), nsmall = 2)
          tempVec[6] <- format(round(median(TempData[TempData$Variable == 0,]$Time), 2), nsmall = 2)
          tempVec[7] <- format(round(median(TempData[TempData$Variable == 1,]$Time), 2), nsmall = 2)
          tempVec[8] <- format(round(mean(TempData[TempData$Variable == 0,]$Time), 2), nsmall = 2)
          tempVec[9] <- format(round(mean(TempData[TempData$Variable == 1,]$Time), 2), nsmall = 2)
          tempVec[10] <- length(unique(TempData[TempData$Variable == 0,]$TMA_ID))
          tempVec[11] <- length(unique(TempData[TempData$Variable == 1,]$TMA_ID))


          medianData <- rbind(medianData, t(tempVec))

          ## Creation of plots
          # Define title for plot
          ifelse(
            !is.null(PlotTitles),
            Title <-
              PlotTitles[i],
            Title <-
              VarName
          )



          # Define colour palette
          PaletteList <-
            c(
              "red",
              "blue",
              "turquoise4",
              "goldenrod4",
              "mediumorchid4",
              "dodgerblue1",
              "green",
              "magenta3"
            )
          PaletteGroups <- nrow(SurvOutputTable)
          Palette <- PaletteList[1:PaletteGroups]


          ## Generate plot
          # Open tiff file
          FileName <- paste(VarName, ".tiff", sep = "")
          tiff(
            FileName,
            # width = 1500,
            # height = 1500,
            width = 2250,
            height = 2500,
            units = "px",
            bg = "transparent",
            res = 300,
            compression = "lzw"
          )

          custom_theme <- function() {
            theme_test() %+replace%
              theme(
                plot.title = element_text(hjust = 0.5, vjust = 1, face = "bold", size = 15)
              )
          }


          suppressMessages(
          plot <- ggsurvplot(
            fit = survfit(Surv(Time, Status) ~ Variable, data = TempData),
            data = TempData,
            linetype = "solid",
            size = 2,
            censor.shape = 124,
            censor.size = 4,
            palette = Palette,
            break.time.by = 12,
            title = Title,
            font.main = c(10, "bold", "Black"),
            xlim = c(0, max(TempData$Time)),
            xlab = LabX,
            font.x = c(10, "bold.italic", "Black"),
            ylab = LabY,
            font.y = c(10, "bold.italic", "Black"),
            legend = "none",
            font.tickslab = 8,
            ggtheme = custom_theme(),
            pval = FALSE,
            pval.size = 5,
            pval.coord = c(0.007, 0.05),
            pval.method = FALSE,
            pval.method.size = 5,
            pval.method.coord = c(0.007, 0.1),
            conf.int = TRUE,
            conf.int.style = "ribbon",
            conf.int.alpha = 0.1,
            # surv.median.line = "hv",
            risk.table = "percentage",
            cumevents = FALSE,
            cumcensor = FALSE,
            tables.col = "strata",
            fontsize = 4,
            risk.table.y.text = FALSE,
            tables.y.text.col = TRUE,
            risk.table.pos = "out"
          ),

          # plot$table <- plot$table +
          #   theme(plot.title = element_text(size = 10, color = "green", face = "bold"))

          plot$table <-
            ggrisktable(
              fit,
              data = TempData,
              color = "strata",
              y.text = TRUE,
              # ylab = "",
              # xlab = "",
              tables.theme = theme_cleantable(),
              font.tickslab = c(20, "bold"),
              fontsize = 2
            )

          )


          # Define Y axis expansion per variable group
          yAxisExp <- 0.2
          VariableLevels <- nrow(SurvOutputTable)
          if (nrow(SurvOutputTable) >= 3) {
            yAxisExp <- 0.2 + (0.15 * (VariableLevels - 2))
          }


          TT1 <- ttheme_minimal(
            core = list(
              bg_params = list(
                fill = Palette,
                alpha = 0.1,
                col = NA
              ),
              fg_params = list(fontface = "bold")
            ),
            colhead = list(fg_params = list(
              col = "Black", fontface = 4L
            )),
            rowhead = list(fg_params = list(
              col = "orange", fontface = 3L
            )),
            base_size = 12
          )

          suppressMessages(
          plot$plot <-
            plot$plot + annotation_custom(
              tableGrob(SurvOutputTable, rows = NULL, theme = TT1),
              xmin = -Inf,
              xmax = Inf,
              ymin = -Inf,
              ymax = 0
            ) +
            scale_y_continuous(expand = c(yAxisExp, 0),
                               breaks = c(0, 0.25, 0.5, 0.75, 1))
          )

          print(plot)

          dev.off()

          zphData <- data.frame(zphP, zphSurv, zphVar)

          zphDataList <- append(zphDataList, list(zphData))


        } else {
          message(paste("The variable '", Variables[i], "' has fewer than two groups."))
        }


        }

        setwd(BSD)
      }


      setwd(OD)

      zphData <- do.call(rbind, zphDataList)
      write.csv(zphData,
                "zphData.csv",)


    }

    colnames(medianData) <-
      c(
        "Cohort",
        "Status",
        "Variable",
        "Median Survival Time (Low)",
        "Median Survival Time (High)",
        "Median of Survival Time (Low)",
        "Median of Survival Time (High)",
        "Mean of Survival Time (Low)",
        "Mean of Survival Time (High)",
        "Low",
        "High"
      )
    csvName <- paste("C:/Program Files/R/Thesis/Data/Survival/medianData", substr(deparse(substitute(Data)), 1, nchar(deparse(substitute(Data))) -5), ".csv", sep = "")
    write.csv(medianData, csvName, row.names = FALSE)

    flexName <- paste("C:/Program Files/R/Thesis/Data/Survival/medianData", substr(deparse(substitute(Data)), 1, nchar(deparse(substitute(Data))) -5), ".docx", sep = "")
    flexMedianData <-
      flextable(medianData) %>%
      add_header_row(colwidths = c(3, 8), values = c("Group Descriptors", "Median Values")) %>%
      align(align = "center", part = "header") %>%
      align(align = "center", part = "body") %>%
      autofit() %>%
      fontsize(i = NULL, j = NULL, size = 7, part = "header") %>%
      fontsize(i = NULL, j = NULL, size = 5.5, part = "body") %>%
      width(width = 0.5) %>%
      width(j = 3, width = 2)

    flexMedianData

    save_as_docx(flexMedianData, path = flexName)

    # Return to base working directory
    setwd(WD)

  }
