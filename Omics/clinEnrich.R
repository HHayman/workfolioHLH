clinEnrich <-
  function(maf = NULL,
           variable = NULL,
           minMut = 5,
           fileName = NULL,
           filePath = NULL) {



    try(enrichData <-
      clinicalEnrichment(maf = maf, clinicalFeature = variable, minMut = minMut))

    enrichData <<- enrichData



    try(plotEnrichmentResults(
      enrich_res = enrichData,
      pVal = 0.05,
      ORthr = 2,
      geneFontSize = 1,
      annoFontSize = 1,
      legendFontSize = 1,
      showTitle = FALSE
    ))



    tiff(
      if (!is.null(fileName)) {
        filename = paste(filePath,
                         fileName,
                         ".tiff",
                         sep = "")
      } else {
        filename = paste(filePath,
                         "clinEnrich",
                         ".tiff",
                         sep = "")
      },
      width = 2000,
      height = 2000,
      units = "px",
      bg = "transparent",
      res = 300,
      compression = "lzw"
    )

    try(plotEnrichmentResults(
      enrich_res = enrichData,
      pVal = 0.05,
      ORthr = 2,
      geneFontSize = 1,
      annoFontSize = 1,
      legendFontSize = 1,
      showTitle = FALSE
    ))

    dev.off()




    setClass("clinEnrich",
             representation(sigGenes = "vector", sigGenesFDR = "vector", sigGenesData = "data.frame", sigGenesFDRData = "data.frame", enrichData = "data.frame"))
    output <-
      new(
        "clinEnrich",
        sigGenes = unique(enrichData$groupwise_comparision[enrichData$groupwise_comparision$p_value < 0.05 &
                                                             enrichData$groupwise_comparision$OR > 2, ]$Hugo_Symbol),
        sigGenesFDR = unique(enrichData$groupwise_comparision[enrichData$groupwise_comparision$p_value < 0.05 &
                                                                enrichData$groupwise_comparision$fdr < 0.05 &
                                                                enrichData$groupwise_comparision$OR > 2,]$Hugo_Symbol),
        sigGenesData = enrichData$groupwise_comparision[enrichData$groupwise_comparision$p_value < 0.05 &
                                                             enrichData$groupwise_comparision$OR > 2, ],
        sigGenesFDRData = enrichData$groupwise_comparision[enrichData$groupwise_comparision$p_value < 0.05 &
                                                                enrichData$groupwise_comparision$fdr < 0.05 &
                                                                enrichData$groupwise_comparision$OR > 2,],
        enrichData = enrichData$groupwise_comparision
      )

    return(output)
  }
