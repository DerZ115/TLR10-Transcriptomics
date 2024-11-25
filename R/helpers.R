packages <- c(
    "DESeq2", "dplyr", "magrittr", "ComplexHeatmap", "tidyverse",
    "RColorBrewer", "stringr", "ggplot2", "ggrepel", "ggpubr",
    "org.Hs.eg.db"
)

lapply(packages, library, character.only = TRUE)

import_data <- function(counts_file, metadata_file, factor_levels) {
    count_data <- read.delim(counts_file) %>%
        # Remove genes with all zero counts
        filter(rowSums(.[, sapply(., is.numeric)]) > 0) %>%
        # Remove version from ENSEMBL id, as it interferes with annotation
        mutate(ENSG = str_remove(ENSG, "\\..*$")) %>%
        # Merge duplicate entries
        group_by(ENSG) %>%
        summarise(across(!where(is.numeric), ~ dplyr::first(.)), across(where(is.numeric), ~ sum(.))) %>%
        ungroup() %>%
        as.data.frame() %>%
        column_to_rownames("ENSG")

    annotation_data <- count_data %>%
        dplyr::select(c(GeneType, GeneSymbol)) %>%
        rownames_to_column("ENSG") %>%
        as_tibble()

    count_data %<>% dplyr::select(-c(GeneType, GeneSymbol))

    sample_data <- read.csv(metadata_file) %>%
        as_tibble() %>%
        mutate(across(names(factor_levels), ~ factor(.x, levels = factor_levels[[cur_column()]]))) %>%
        mutate(group = factor(paste(celltype, condition, time, sep = "_")))

    # Check if Sample names in countdata and sampledata match
    # otherwise reorder countdata accordingly
    if (!all(colnames(count_data) == sample_data$sampleID)) {
        countdata <- count_data[, sample_data$sampleID]
    }

    list(count_data = count_data, annotation_data = annotation_data, sample_data = sample_data)
}
