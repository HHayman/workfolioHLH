groupLadder <-
  function(mafA = NULL,
           mafB = NULL,
           mafNameA = NULL,
           mafNameB = NULL,
           genes = NULL,
           ID = NULL,
           ySize = NULL,
           pathname = NULL,
           outputname = NULL) {
    #### Configs ####
    allID <-
      unique(c(
        mafA@data$Tumor_Sample_Barcode,
        mafB@data$Tumor_Sample_Barcode
      ))
    nameID <- ID

    template <-
      function(IDs = NULL,
               gene = NULL,
               mafName = NULL) {
        templateData <-
          data.frame(
            ID = IDs,
            Gene = gene,
            TotalMut = "Widltype",
            maf = mafName
          )
      }

    `%!in%` <- Negate(`%in%`)



    #### Parse data ####
    mutDataA <- list()
    for (i in 1:length(genes)) {
      mutDataTemp <- geneStatus(maf = mafA,
                                ID = ID,
                                genes = genes[i]) %>% select(ID, Gene, TotalMut) %>%
        mutate(maf = mafNameA, TotalMut = "Mutant")
      mutDataTemp <- rbind(
        mutDataTemp,
        template(
          IDs = allID,
          gene = genes[i],
          mafName = mafNameA
        ) %>% filter(ID %!in% mutDataTemp$ID)
      )
      mutDataA[[i]] <- mutDataTemp
    }

    mutDataB <- list()
    for (i in 1:length(genes)) {
      mutDataTemp <- geneStatus(maf = mafB,
                                ID = ID,
                                genes = genes[i]) %>% select(ID, Gene, TotalMut) %>%
        mutate(maf = mafNameB, TotalMut = "Mutant")
      mutDataTemp <- rbind(
        mutDataTemp,
        template(
          IDs = allID,
          gene = genes[i],
          mafName = mafNameB
        ) %>% filter(ID %!in% mutDataTemp$ID)
      )
      mutDataB[[i]] <- mutDataTemp
    }

    mutData <- do.call(rbind, c(mutDataA, mutDataB)) %>%
      dplyr::rename(Status = TotalMut) %>%
      mutate(Group = paste(Gene, "\n", maf, sep = "")) %>%
      select(Status, ID, Group) %>%
      arrange(ID, Group) %>%
      dplyr::mutate(Paired = rep(1:(n() / 2), each = 2),
                    Group = factor(Group))

    tempIDsGene <- list()
    for (i in 1:length(genes)) {
      mutDataTemp <-
        mutData %>% filter(grepl(genes[i], Group, fixed = TRUE))
      tempIDs <- unique(mutDataTemp$ID)
      tempIDsList <- list()
      for (j in 1:length(unique(mutDataTemp$ID))) {
        tempChange <- mutDataTemp %>%
          filter(ID == unique(mutDataTemp$ID)[j]) %>%
          mutate(fieldChange = ifelse(length(unique(Status)) == 1, "No Change", "Change")) %>%
          select(ID, Group, fieldChange)
        tempIDsList[[j]] <- tempChange
      }
      tempIDsGene[[i]] <- do.call(rbind, tempIDsList)
    }
    fcKey <- do.call(rbind, tempIDsGene)
    mutData <- full_join(mutData, fcKey, by = c("ID", "Group")) %>%
      mutate(fieldChangeLine = ifelse(fieldChange == "No Change", "solid", "dotted"))



    #### Generate plots ####
    tiff(
      filename = paste(pathname,
                       ifelse(is.null(outputname), "groupLadder.tiff", outputname),
                       sep = ""),
      width = 4000,
      height = 4000,
      units = "px",
      bg = "transparent",
      res = 300,
      compression = "lzw"
    )
    print(
      ggplot(mutData, aes(x = Group,
                          y = ID)) +
        geom_line(aes(group = Paired, linetype = fieldChange)) +
        geom_point(aes(color = Status), size = 4) +
        labs(x = "Group", y = nameID) +
        theme_classic() +
        theme(
          legend.position = "right",
          legend.text = element_text(size = 14, face = "bold"),
          legend.title = element_text(size = 18, face = "bold"),
          axis.title.x = element_text(size = 18, face = "bold"),
          axis.text.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 18, face = "bold"),
          axis.text.y = element_text(size = ifelse(!is.null(ySize), ySize, 0), face = "bold")
        ) +
        scale_linetype_manual(
          name = "Field Change",
          values = c("solid", "dotted"),
          labels = c("Change", "No Change")
        )
    )
    dev.off()



    return(mutData)


  }
