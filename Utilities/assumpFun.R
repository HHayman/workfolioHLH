
## Packages
library(dplyr)
library(ggpubr)
library(tidyr)



assumpFun <-
  function(localData = NULL,
           columns = NULL,
           pathName = NULL) {
    for (i in 1:length(columns)) {
      ## qqPlots
      tiff(
        filename = paste(
          pathName,
          "qqPlots/",
          columns[i],
          ".tiff",
          sep = ""
        ),
        width = 4000,
        height = 4000,
        units = "px",
        bg = "transparent",
        res = 300,
        compression = "lzw"
      )
      print(ggqqplot(localData, x = columns[i]))
      dev.off()

      ## Histograms
      tiff(
        filename = paste(
          pathName,
          "Histograms/",
          columns[i],
          ".tiff",
          sep = ""
        ),
        width = 4000,
        height = 4000,
        units = "px",
        bg = "transparent",
        res = 300,
        compression = "lzw"
      )
      print(ggplot(localData, aes(x = .data[[columns[i]]])) + geom_histogram() + theme_classic())
      dev.off()
    }

    ## Normality tests
    normalityTests <- localData %>%
      select(columns) %>%
      pivot_longer(columns, names_to = "Marker", values_to = "Expression") %>%
      group_by(Marker) %>%
      summarise_all(.funs = funs(
        statistic = shapiro.test(.)$statistic,
        p.value = shapiro.test(.)$p.value
      )) %>%
      mutate(
        sig = ifelse(p.value < 0.05, "Sig", "Not sig"),
        normal = ifelse(p.value < 0.05, "Not normal", "Normal")
      )
    write.csv(
      normalityTests,
      paste(pathName, "shapiroTests.csv", sep = ""),
      row.names = FALSE
    )



  }
