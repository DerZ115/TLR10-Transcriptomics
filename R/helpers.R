packages <- c(
    "DESeq2", "dplyr", "magrittr", "ComplexHeatmap", "tidyverse",
    "RColorBrewer", "stringr", "ggplot2", "ggrepel", "ggpubr", "ggforce",
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
        mutate(class = factor(paste(celltype, condition, time, sep = "_")))

    # Check if Sample names in countdata and sampledata match
    # otherwise reorder countdata accordingly
    if (!all(colnames(count_data) == sample_data$sampleID)) {
        countdata <- count_data[, sample_data$sampleID]
    }

    list(count_data = count_data, annotation_data = annotation_data, sample_data = sample_data)
}

plot_pca <- function(data, color_var, shape_var, n = 500, maxPC = 2, center = TRUE, scale = TRUE, 
                     plot_all = F) {
  samples <- colData(data) %>% 
    as.data.frame()
  
  
  data.pca <- assay(data) %>%
    t() %>%
    prcomp(center = center, scale. = scale) 
  
  data.pca.plot <- data.pca$x %>%
    as.data.frame() %>%
    rownames_to_column("sampleID") %>%
    left_join(samples, by = "sampleID")
  
  data.pca.var <- as.data.frame(t(summary(data.pca)$importance)) %>%
    set_colnames(c("sd", "var_prop", "var_cumprop")) %>%
    mutate(across(c(var_prop, var_cumprop), ~ .x * 100), Component = seq.int(1, nrow(.), 1))
  
  p1 <- ggplot(data.pca.var, aes(x = Component)) +
    geom_bar(aes(y = var_prop), stat = "identity") +
    geom_line(aes(y = var_cumprop)) +
    geom_point(aes(y = var_cumprop)) +
    scale_x_continuous(limits = c(0.5, nrow(data.pca.var) + 0.5), 
                       breaks = seq(1, nrow(data.pca.var), by = 1),
                       expand = expansion(add=0)) +
    scale_y_continuous(expand = expansion(mult=c(0, 0.05))) +
    ylab("Variance (%)") +
    theme_pubr() +
    theme(panel.grid.major.y = element_line(linewidth = 0.5))
  
  print(p1)
  
  
  p2 <- ggplot(data.pca.plot, aes(x = PC1, y = PC2, 
                                  color = interaction(!!!syms(color_var), sep=":", lex.order = T),
                                  shape = interaction(!!!syms(shape_var), sep=":", lex.order = T))) +
    geom_point(size = 2) +
    labs(color = "") +
    theme_pubr() +
    theme(panel.grid.major = element_line(linewidth = 0.5),
          legend.position = "right")
  
  print(p2)
  
  if (plot_all) {
    
    p3 <- ggplot(data.pca.plot, aes(color = interaction(!!!syms(color_var), sep=":", lex.order = T),
                                    shape = interaction(!!!syms(shape_var), sep=":", lex.order = T))) +
      geom_autopoint(size = 2) +
      geom_autodensity(position = "identity", fill = NA) +
      facet_matrix(vars(paste0("PC", seq.int(1, maxPC, 1))), layer.diag = 2) +
      labs(color = "") +
      theme_bw() +
      theme(panel.grid.major = element_line(linewidth = 0.5))
    
    print(p3)
  }
}

plot_counts <- function(dds, annotation_data, counts_to_plot, x_var, color_var, genename_col = "GeneSymbol", 
id_col = "ENSG") {
    counts.plot <- lapply(counts_to_plot, function(x) {
      plotCounts(dds, annotation_data[(annotation_data %>% pull(genename_col)) == x,] %>% 
                 pull(id_col), intgroup = c(x_var, color_var), returnData = TRUE)
  }) %>%
  set_names(counts_to_plot) %>%
  lapply(rownames_to_column, var = "Sample_ID") %>%
  bind_rows(.id = genename_col) %>%
  group_by(across(all_of(c(genename_col, x_var, color_var)))) %>%
  summarise(mean_count = mean(count), sd_count = sd(count), .groups = "drop") %>%
  mutate(min_error = ifelse(mean_count - sd_count >= 0, mean_count - sd_count, 0), 
           max_error = mean_count + sd_count)

ggplot(counts.plot, aes(x = interaction(!!!syms(x_var), sep=":", lex.order = T), 
                        y = mean_count, fill = interaction(!!!syms(color_var), sep=":", lex.order = T))) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = min_error, ymax = max_error), 
                position = position_dodge(0.9), width = 0.2) +
  labs(x = "", y = "Mean Count", fill = "", title = "TLR10 Expression") +
  theme_pubr() +
  theme(legend.position = "right", legend.text = element_text(size = 12), panel.grid.major.y = element_line(linewidth = 0.5))
}

plot_heatmap <- function(data, color_scale, column_colors, genedata, samples, 
                         pvalues, title = "Rel. Expression", cluster_rows = TRUE, 
                         cluster_cols = TRUE, max_rows = 20, scale_by_row = TRUE) {
  
  annotation.row <- genedata %>%
    filter(ENSG %in% rownames(data)) 
  
  annotation.col <- samples %>%
    filter(sampleID %in% colnames(data)) %>%
    dplyr::select(c(class))
  
  annotation_labels <- unique(annotation.col$class)
  annotation_colors <- column_colors[annotation_labels]
  
  if (scale_by_row) {
    data %<>%
      t() %>%
      scale() %>%
      t()
  }
  
  if (nrow(data) <= max_rows) {
    hm <- Heatmap(as.matrix(data), col = color_scale,  
                  name = title, cluster_rows = cluster_rows, column_split = annotation.col$class,
                  show_row_dend = TRUE, show_column_names = FALSE, row_labels = annotation.row$GeneSymbol,
                  show_row_names = TRUE, row_title = NULL, column_title = NULL, cluster_columns = cluster_cols, 
                  top_annotation = HeatmapAnnotation(class = anno_block(labels = annotation_labels,
                                                                        gp = gpar(fill = annotation_colors))),
                  heatmap_legend_param = list(legend_height = unit(10, "cm"),
                                              legend_width = unit(1.5, "cm"),
                                              title_position = "lefttop-rot"))
    
    draw(hm)
    
  } else {
    hm <- Heatmap(as.matrix(data), col = color_scale,  
                  name = title, cluster_rows = cluster_rows, column_split = annotation.col$class, 
                  show_row_dend = FALSE, show_column_names = FALSE, show_row_names = FALSE, 
                  row_title = NULL, column_title = NULL, cluster_columns = cluster_cols, 
                  top_annotation = HeatmapAnnotation(class = anno_block(labels = annotation_labels,
                                                                        gp = gpar(fill = annotation_colors))),
                  heatmap_legend_param = list(legend_height = unit(10, "cm"),
                                              legend_width = unit(1.5, "cm"),
                                              title_position = "lefttop-rot"))
    
    draw(hm)
    
    data.top <- data %>% 
      as.data.frame() %>%
      rownames_to_column("ENSG") %>%
      left_join(data.frame(pvalues, ENSG = names(pvalues)), by = "ENSG") %>%
      column_to_rownames("ENSG") %>%
      arrange(pvalues) %>%
      dplyr::select(!pvalues) %>%
      head(max_rows) 
    
    annotation.row.top <- annotation.row %>%
      filter(ENSG %in% rownames(data.top))
    
    hm <- Heatmap(as.matrix(data.top), col = color_scale,  
                  name = title, cluster_rows = cluster_rows, column_split = annotation.col$class, 
                  show_row_dend = TRUE, show_column_names = FALSE, row_labels = annotation.row.top$Name,
                  show_row_names = TRUE, row_title = NULL, column_title = NULL, cluster_columns = cluster_cols, 
                  top_annotation = HeatmapAnnotation(class = anno_block(labels = annotation_labels,
                                                                        gp = gpar(fill = annotation_colors))),
                  heatmap_legend_param = list(legend_height = unit(10, "cm"),
                                              legend_width = unit(1.5, "cm"),
                                              title_position = "lefttop-rot"))
    
    draw(hm)
    
  }
}
