#' #'CutFlow
#' 
#' 
#' #' @author Hannah Hayman, \email{hannah.louise.hayman@@gmail.com}
#' #' This function takes in a dataset(s), produces cut points using the maximally selected rank statistic and recodes your dataset(s).
#' 
#' 
#' #' @import ggplot2
#' #' @import survminer
#' #' @import gplots
#' 
#' 
#' #' @param Subdirectory A subdirectory of your working directory, in which you store your dataset CSV files. Required.
#' #' @param TrainingData The name of your training dataset, which is used to generate cutpoints. Required.
#' #' @param CutPointStatus The survival status you wish to use to define your cutpoints. Required.
#' #' @param CutPointTime The time variable you wish to use to define your cutpoints. Required.
#' #' @param minprop The minimum proportion of cases to be included in a group. Default is 0.1 if none is provided.
#' #' @param Variables A list of all variables you wish to be included. Names must match and be identical across datasets. Required.
#' #' @param Greyscale Use 'TRUE' if a greyscale plot is wanted. Otherwise, red/blue is the default palette.
#' 
#' 
#' #' @return CutFlow returns a folder with subdirectories, containing your input datasets, the output of cutpoint generation and coded versions of your datasets.
#' 
#' 
#' #' @examples CutFlow(Subdirectory = "CutFlowData", TrainingData = "Glasgow", CutPointStatus = "CSS", CutPointTime = "CSS_2017", minprop = 0.1, Greyscale = TRUE, Variables = c( "GD_PercPositiveCellsinHealthyTissue", "GD_PercPositiveCellsinTumourTissue",  "CD8_PercPositiveCellsinHealthyTissue", "CD8_PercPositiveCellsinTumourTissue"))
#' #' @export




# Packages
library(ggplot2)
library(survminer)
library(gplots)



# Function
CutFlow <-
  function(Subdirectory,
           TrainingData,
           CutPointStatus,
           CutPointTime,
           minprop = 0.1,
           Variables,
           Greyscale = FALSE) {
    # Read data and check for errors
    #---------------------------------------------------------------------------



    # Save original working directory and set to restore on exit
    WD <- toString(getwd())
    on.exit(setwd(WD), add = TRUE)

    # Check that the given subdirectory for reading data exists. Error if not
    if (file.exists(Subdirectory)) {
      setwd(Subdirectory)
    } else {
      stop(
        'Incompatible subdirectory. Please check the name of the folder containing your files.It should be a subdirectory within your working directory.'
      )
    }

    # Read and name data files
    filenames <-
      gsub("\\.csv$", "", list.files(pattern = "\\.csv$"))
    Datasets <-
      lapply(
        list.files(pattern = "\\.csv$"),
        read.csv,
        sep = ",",
        fileEncoding = 'UTF-8-BOM'
      )
    names(Datasets) <- filenames
    
    # Revert to original working directory
    setwd(WD)

    # Check that datasets are dataframes
    for (i in 1:length(filenames)) {
      if (!inherits(Datasets[[filenames[i]]], "data.frame"))
        stop("data should be an object of class data.frame")
      Datasets[[filenames[i]]] <-
        as.data.frame(Datasets[[filenames[i]]])
    }

    # Subset the training dataset
    TrainingDataSetName <- toString(TrainingData)
    TrainingDataSet <- Datasets[[TrainingDataSetName]]

    # Check that cut point variables are in the training dataset
    if (!all(c(CutPointStatus, CutPointTime) %in% colnames(TrainingDataSet)))
      stop("Specify correct column names containing CutPointStatus and CutPointTime values.")




    # Create a subdirectory for output
    #---------------------------------------------------------------------------

    # Define and check subdirectory name
    Number <- 1
    OutputDirectory <-
      paste("CutFlow", Sys.Date(), Number, sep = '_')
    while (file.exists(OutputDirectory))
    {
      Number <- Number + 1
      OutputDirectory <-
        paste("CutFlow", Sys.Date(), Number, sep = '_')
    }

    # Create the subdirectory
    dir.create(OutputDirectory)
    setwd(OutputDirectory)
    OD <- toString(getwd())






    # Deposit raw datasets into a folder for record keeping
    #---------------------------------------------------------------------------

    # Create and set a subdirectory
    dir.create("0.OriginalDatasetFiles")
    setwd("0.OriginalDatasetFiles")

    # Write datasets to csv files
    for (i in 1:length(filenames)) {
      write.csv(Datasets[[filenames[i]]], paste(filenames[i], ".csv", sep = ""), row.names = FALSE)
    }

    #Revert to output subdirectory
    setwd(OD)







    # Define cutpoints using the 'survminer' package
    #---------------------------------------------------------------------------
    # Create and set a subdirectory
    dir.create("1.CutPointOutput")
    setwd("1.CutPointOutput")

    # Record user inputs
    TrainingData_Record <- TrainingData
    CutPointStatus_Record <- CutPointStatus
    CutPointTime_Record <- CutPointTime
    MinProp_Record <- minprop
    Variables_Record <- Variables

    Syntax_Record <-
      data.frame(
        TrainingData_Record,
        CutPointStatus_Record,
        CutPointTime_Record,
        MinProp_Record,
        Variables_Record
      )

    write.csv(Syntax_Record,
              paste("Syntax_Record.csv", sep = ""),
              row.names = FALSE)


    # Open a sink for console messages
    sink("waste.txt")

    # Create empty lists to hold output
    Plots <- list()
    CutPoints <- list()

    # Define the number of variables to produce a cutpoint for
    n = length(Variables)


    ## Generate and plot the cutpoints

    # Open a loop to cycle through variables
    for (i in 1:n) {
      # Define filenames for cutpoint graphs
      TiffName <-
        paste(TrainingDataSetName,
              "_",
              CutPointStatus,
              "_",
              Variables[i],
              ".tiff",
              sep = "")
      # Define an empty TIFF file for storing output
      tiff(TiffName)

      # Generate cutpoints
      NewCut <-
        surv_cutpoint(
          Datasets[[TrainingData]],
          time = CutPointTime,
          event = CutPointStatus,
          minprop = minprop,
          Variables[i]
        )

      # Set colour palettes
      ifelse (
        Greyscale == TRUE,
        PlotPalette <-
          c("black", "slategrey"),
        PlotPalette <- c("#d70033", "#5596e6")
      )

      # Plot the newly generated cutpoints
      NewPlot <-
        plot(NewCut,
             Variables[i],
             palette = c(PlotPalette),
             main = CutPointStatus)

      # Print the cutpoint plot to the empty TIFF file
      print(NewPlot)
      dev.off()

      # Add plots to plots list
      Plots <- c(Plots, NewPlot)

    }

    # Generate cutpoints
    CutPoints <-
      surv_cutpoint(
        Datasets[[TrainingData]],
        time = CutPointTime,
        event = CutPointStatus,
        minprop = minprop,
        Variables
      )

    # Summarise cutpoints
    CutPoints <- summary(CutPoints)

    # Create a vector of cutpoints
    CutPointVector <- as.vector(CutPoints[["cutpoint"]])

    # Define title for cutpoints file
    CutPointsTitle <-
      paste(TrainingDataSetName,
            " ",
            CutPointStatus,
            " ",
            Variables[i],
            " Cutpoints",
            sep = "")

    # Define name for cutpoints file
    PDFname <-
      paste(TrainingDataSetName, "_", CutPointStatus, ".pdf", sep = "")

    # Define an empty PDF file for storing output
    pdf(PDFname)

    # Print cutpoints to PDF
    textplot(CutPoints,
             halign = "center",
             valign = "top",
             cex = 0.6)
    # Add title to PDF
    title(paste(TrainingDataSetName,
                "_",
                CutPointStatus,
                "_Cutpoints",
                sep = ""))

    # Print the plots to the PDF and close the PDF
    print(Plots)
    dev.off()

    # Close the sink for console messages
    sink()

    #Revert to output subdirectory
    setwd(OD)





    # Apply cutpoints to the datasets
    #---------------------------------------------------------------------------

    # Create and set a subdirectory
    dir.create("2.CodedDatasets")
    setwd("2.CodedDatasets")

    # Create copy of datasets for coding
    CodedDatasets <- Datasets

    # Add empty variable to all datasets, suffix '_coded'
    for (k in 1:length(CodedDatasets)) {
      for (i in 1:length(Variables)) {
        if (Variables[i] %in% colnames(CodedDatasets[[k]])) {
          CodingVariable <- paste(Variables[i], "_Coded", sep = "")
          DS <- CodedDatasets[[k]]
          Var <- DS[[Variables[i]]]
          CP <- CutPoints[i, 1]
          CodedDatasets[[k]][paste(CodingVariable)] <-
            ifelse(Var <= CP, 0, 1)
        }
      }
      CodedDatasets[[k]][is.na(CodedDatasets[[k]])] <- ""
    }

    # Write datasets to csv files
    CodedNames <- names(CodedDatasets)
    for (i in 1:length(CodedNames)) {
      write.csv(CodedDatasets[[CodedNames[i]]],
                paste(CodedNames[i], ".csv", sep = ""),
                row.names = FALSE)


    }

    # Revert to working directory
    setwd(WD)

  }
