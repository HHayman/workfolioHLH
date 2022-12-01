library(ggplot2)
library(ggrepel)


intSurvLymph <-
  function(lymphData = NULL,
           survData = NULL,
           cellType = NULL,
           interactions = "No",
           reference = NULL,
           labelSize = 3,
           titleSize = 12,
           survVar = "Prognosis",
           filePath = NULL) {

    if(is.null(reference)){reference <- "NA"}

    #### Lymph Plot ####
    lymphData$Significance <-
      ifelse(
        lymphData$p_value <= 0.001,
        "< 0.001",
        ifelse(
          lymphData$p_value < 0.01,
          "< 0.01",
          ifelse(lymphData$p_value < 0.05, "< 0.05", "NS")
        )
      )
    lymphData$Direction <-
      ifelse(lymphData$OR >= 2,
             "Increased",
             ifelse(lymphData$OR <= 0.5, "Decreased", "No Difference"))

    lymphData$p_value <- ifelse(lymphData$p_value == 0.000, 0.001, lymphData$p_value)
    lymphData$OR <- ifelse(lymphData$OR == 0, min(lymphData[lymphData$OR != 0,]$OR), lymphData$OR)
    lymphData$OR <- ifelse(lymphData$OR == "Inf", max(lymphData[lymphData$OR != "Inf",]$OR), lymphData$OR)

    if(interactions == "Yes"){
    lymphData <- lymphData[str_count(lymphData$Hugo_Symbol, "_") == 1, ]
    }

    checkLymph <<- lymphData

    lymphPlot <-
      ggplot(lymphData,
             aes(
               x = log10(OR),
               y = -log10(p_value),
               colour = Direction,
               size = Significance
             )) +
      xlab("log(Odds Ratio)") +
      ylab("-log(P Value)") +
      geom_point() +
      theme_classic() +
      theme(plot.title = element_text(
        face = "bold",
        colour = "black",
        size = titleSize,
        hjust = 0.5
      ), legend.position = "right", legend.text = element_text(size = 12)) +
      geom_hline(
        yintercept = -log10(0.05),
        linetype = "dashed",
        colour = "grey",
        size = 0.5
      ) +
      geom_vline(
        xintercept = log10(2),
        linetype = "dashed",
        colour = "grey",
        size = 0.5
      ) + annotate(
        x = log10(2),
        y = max(-log10(lymphData$p_value)),
        label = "OR = 2",
        geom = "label",
        size = 3
      ) +
      geom_vline(
        xintercept = log10(0.5),
        linetype = "dashed",
        colour = "grey",
        size = 0.5
      ) + annotate(
        x = log10(0.5),
        y = max(-log10(lymphData$p_value)),
        label = "OR = 0.5",
        geom = "label",
        size = 3
      ) +
      annotate(
        y = -log10(0.05),
        x = max(log10(lymphData$OR)),
        label = "P Value \n= 0.05",
        geom = "label",
        size = 3
      ) +
      geom_text_repel(
        data = subset(lymphData, p_value < 0.05 &
                        OR >= 2 | p_value < 0.05 & OR <= 0.5),
        aes(
          x = log10(OR),
          y = -log10(p_value),
          label = Hugo_Symbol
        ),
        show.legend = FALSE,
        size = labelSize
      ) +
      scale_color_manual(
        name = paste("Likelihood of Mutation\n(Ref: ", reference, ")", sep = ""),
        breaks = c("Increased", "Decreased", "No Difference"),
        values = c("red", "deepskyblue", "black")
      ) +
      scale_size_manual(
        name = "Significance",
        breaks = c("NS", "< 0.05", "< 0.01", "< 0.001"),
        values = c(1, 2, 3, 4)
      )






    #### Surv Plot ####
    survData$P_value <- as.numeric(survData$P_value)
    survData$HR <- as.numeric(survData$HR)
    survData$Significance <-
      ifelse(
        survData$P_value <= 0.001,
        "< 0.001",
        ifelse(
          survData$P_value < 0.01,
          "< 0.01",
          ifelse(survData$P_value < 0.05, "< 0.05", "NS")
        )
      )
    survData$Direction <-
      ifelse(survData$HR >= 2,
             "Unfavourable",
             ifelse(survData$HR <= 0.5, "Favourable", NA))

    survData$P_value <- ifelse(survData$P_value == 0.000, 0.001, survData$P_value)
    survData$HR <- ifelse(survData$HR == 0, min(survData[survData$HR != 0,]$HR), survData$HR)

    checkSurv1 <<- survData

    if(interactions == "Yes"){
    survData <- survData[grep("_", survData$Genes_combinations),]
    }

    checkSurv2 <<- survData

    survPlot <-
      ggplot(survData,
             aes(
               x = log10(HR),
               y = -log10(P_value),
               colour = Direction,
               size = Significance
             )) +
      xlab("log(Hazard Ratio)") +
      ylab("-log(P Value)") +
      geom_point() +
      theme_classic() +
      theme(plot.title = element_text(
        face = "bold",
        colour = "black",
        size = titleSize,
        hjust = 0.5
      ), legend.position = "right", legend.text = element_text(size = 12)) +
      geom_hline(
        yintercept = -log10(0.05),
        linetype = "dashed",
        colour = "grey",
        size = 0.5
      ) +
      geom_vline(
        xintercept = log10(2),
        linetype = "dashed",
        colour = "grey",
        size = 0.5
      ) + annotate(
        x = log10(2),
        y = max(-log10(survData[survData$HR != "Inf",]$P_value)),
        label = "HR = 2",
        geom = "label",
        size = 3
      ) +
      geom_vline(
        xintercept = log10(0.5),
        linetype = "dashed",
        colour = "grey",
        size = 0.5
      ) + annotate(
        x = log10(0.5),
        y = max(-log10(survData$P_value)),
        label = "HR = 0.5",
        geom = "label",
        size = 3
      ) +
      annotate(
        y = -log10(0.05),
        x = max(log10(survData$HR)),
        label = "P Value \n= 0.05",
        geom = "label",
        size = 3
      ) +
      geom_text_repel(
        data = subset(survData, P_value < 0.05 &
                        HR >= 2 | P_value < 0.05 & HR <= 0.5),
        aes(
          x = log10(HR),
          y = -log10(P_value),
          label = Genes_combinations
        ),
        show.legend = FALSE,
        size = labelSize
      ) +
      scale_colour_manual(
        name = "Prognosis",
        breaks = c("Unfavourable", "Favourable"),
        values = c("red", "deepskyblue")
      ) +
      scale_size_manual(
        name = "Significance",
        breaks = c("NS", "< 0.05", "< 0.01", "< 0.001"),
        values = c(1, 2, 3, 4)
      )







    ### Dual Plot ####
    survData <-
      dplyr::rename(survData,
                    survP = P_value,
                    survHR = HR,
                    Hugo_Symbol = Genes_combinations)
    survData <- survData[c("Hugo_Symbol", "survP", "survHR")]
    lymphData <-
      dplyr::rename(lymphData, lymphP = p_value, lymphOR = OR)
    lymphData <- lymphData[c("Hugo_Symbol", "lymphP", "lymphOR")]
    dualData <- merge(survData, lymphData, by = "Hugo_Symbol")

    checkDual1 <<- dualData

    dualData$codeSurvHR <-
      ifelse(dualData$survHR >= 2,
             "Unfavourable",
             ifelse(dualData$survHR <= 0.5, "Favourable", NA))

    dualData$codeLymphOR <-
      ifelse(dualData$lymphOR >= 2,
             "Increased",
             ifelse(dualData$lymphOR <= 0.5, "Decreased", "No Difference"))

    dualData$SignificanceS <-
      ifelse(
        dualData$survP <= 0.001,
        "< 0.001",
        ifelse(
          dualData$survP < 0.01,
          "< 0.01",
          ifelse(dualData$survP < 0.05, "< 0.05", "NS")
        )
      )

    checkDual2 <<- dualData

    dualData$SignificanceL <-
      ifelse(
        dualData$lymphP <= 0.001,
        "< 0.001",
        ifelse(
          dualData$lymphP < 0.01,
          "< 0.01",
          ifelse(dualData$lymphP < 0.05, "< 0.05", "NS")
        )
      )

    dualData$Significance <- ifelse(dualData$SignificanceS != "NS" & dualData$SignificanceL != "NS", "Significant", "Not Significant")

    if(interactions == "Yes"){
    dualData <- dualData[str_count(dualData$Hugo_Symbol, "_") == 1, ]
}


    checkDual <<- dualData


    dualPlot <-
      ggplot(dualData,
             aes(
               x = log10(lymphOR),
               y = log10(survHR),
               colour = codeSurvHR,
               size = Significance,
               shape = codeLymphOR
             )) +
      xlab("log(Odds Ratio)") +
      ylab("log(Hazard Ratio)") +
      geom_point() +
      theme_classic() +
      theme(plot.title = element_text(
        face = "bold",
        colour = "black",
        size = titleSize,
        hjust = 0.5
      ), legend.position = "right", legend.text = element_text(size = 12)) +
      geom_vline(
        xintercept = log10(2),
        linetype = "dashed",
        colour = "grey",
        size = 0.5
      ) +
        annotate(
        x = log10(2),
        y = max(log10(dualData$survHR)),
        label = "OR = 2",
        geom = "label",
        size = 3
      ) +
      geom_vline(
        xintercept = log10(0.5),
        linetype = "dashed",
        colour = "grey",
        size = 0.5
      ) +
      annotate(
        x = log10(0.5),
        y = max(log10(dualData$survHR)),
        label = "OR = 0.5",
        geom = "label",
        size = 3
      ) +
      geom_hline(
        yintercept = log10(2),
        linetype = "dashed",
        colour = "grey",
        size = 0.5
      ) +
      annotate(
        y = log10(2),
        x = max(log10(dualData$lymphOR)),
        label = "OR = 2",
        geom = "label",
        size = 3
      ) +
      geom_hline(
        yintercept = log10(0.5),
        linetype = "dashed",
        colour = "grey",
        size = 0.5
      ) +
      annotate(
        y = log10(0.5),
        x = max(log10(dualData$lymphOR)),
        label = "OR = 0.5",
        geom = "label",
        size = 3
      ) +
      geom_text_repel(
        data = subset(dualData, survP < 0.05 & survHR >= 2 | survP < 0.05 & survHR <= 0.5),
        aes(
          x = log10(lymphOR),
          y = log10(survHR),
          label = Hugo_Symbol
        ),
        show.legend = FALSE,
        size = labelSize
      ) +
      scale_colour_manual(
        name = "Prognosis",
        breaks = c("Unfavourable", "Favourable"),
        values = c("red", "deepskyblue")
      ) +
      scale_size_manual(
        name = "Significance",
        breaks = c("Not Significant", "Significant"),
        values = c(1, 2)
      ) +
      scale_shape_manual(
        name = paste("Likelihood of Mutation\n(Ref: ", reference, ")", sep = ""),
        breaks = c("Increased", "Decreased", "No Difference"),
        values = c(24, 25, 16)
      )





    #### Outputs ####
    ## Lymph plot
    tiff(
      filename = paste(filePath,
                       "lymphPlot.tiff",
                       sep = ""),
      width = 2500,
      height = 2000,
      units = "px",
      bg = "transparent",
      res = 300,
      compression = "lzw"
    )
    print(lymphPlot)
    dev.off()
    print(lymphPlot)



    ## Surv plot
    tiff(
      filename = paste(filePath,
                       "survPlot.tiff",
                       sep = ""),
      width = 2500,
      height = 2000,
      units = "px",
      bg = "transparent",
      res = 300,
      compression = "lzw"
    )
    print(survPlot)
    dev.off()
    print(survPlot)



    ## Dualplot
    tiff(
      filename = paste(filePath,
                       "dualPlot.tiff",
                       sep = ""),
      width = 2500,
      height = 2000,
      units = "px",
      bg = "transparent",
      res = 300,
      compression = "lzw"
    )
    print(dualPlot)
    dev.off()
    print(dualPlot)




    return(list(lymphPlot = lymphPlot, survPlot = survPlot, dualPlot = dualPlot))



  }
