packages <- c(
    "DESeq2", "dplyr", "magrittr", "ComplexHeatmap", "tidyverse",
    "RColorBrewer", "stringr", "ggplot2", "ggrepel", "ggpubr", "ggforce",
    "org.Hs.eg.db"
)

lapply(packages, library, character.only = TRUE)

#' Import Data
#'
#' This function imports count data and metadata, processes the count data by removing genes with all zero counts,
#' removing version numbers from ENSEMBL IDs, merging duplicate entries, and reordering columns to match sample data.
#' It also processes the metadata by converting specified columns to factors with given levels and creating a new
#' 'class' column by combining cell type, condition, and time.
#'
#' @param counts_file Path to the counts file.
#' @param metadata_file Path to the metadata file.
#' @param factor_levels List of factor levels for the metadata.
#' @return A list containing count data, annotation data, and sample data.
#' @export
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

#' Plot PCA
#'
#' This function performs Principal Component Analysis (PCA) on the DESeq2 dataset and generates PCA plots.
#' It creates a bar plot showing the variance explained by each principal component and a scatter plot of the first
#' two principal components, colored and shaped by specified variables. Optionally, it can plot all principal components.
#'
#' @param data DESeq2 dataset.
#' @param color_var Variables to color by.
#' @param shape_var Variables to shape by.
#' @param n Number of top genes to use.
#' @param maxPC Maximum number of principal components to plot.
#' @param center Whether to center the data.
#' @param scale Whether to scale the data.
#' @param plot_all Whether to plot all principal components.
#' @return None.
#' @export
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

#' Plot Counts
#'
#' This function generates a bar plot of mean counts for specified genes, grouped by specified variables.
#' It calculates mean and standard deviation of counts, and adds error bars to the plot.
#'
#' @param dds DESeq2 dataset.
#' @param annotation_data Annotation data.
#' @param counts_to_plot Genes to plot.
#' @param x_var Variables for x-axis.
#' @param color_var Variables for color.
#' @param genename_col Column name for gene names.
#' @param id_col Column name for gene IDs.
#' @return A ggplot object.
#' @export
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

#' Plot Heatmap
#'
#' This function generates a heatmap of expression data, with options to cluster rows and columns, scale data by row,
#' and limit the number of rows plotted. It also adds annotations for columns and rows based on provided gene and sample data.
#'
#' @param data Expression data.
#' @param color_scale Color scale for the heatmap.
#' @param column_colors Colors for columns.
#' @param column_labels Labels for columns.
#' @param genedata Gene data.
#' @param samples Sample data.
#' @param pvalues P-values for the genes.
#' @param title Title of the heatmap.
#' @param cluster_rows Whether to cluster rows.
#' @param cluster_cols Whether to cluster columns.
#' @param max_rows Maximum number of rows to plot.
#' @param scale_by_row Whether to scale by row.
#' @return None.
#' @export
plot_heatmap <- function(data, color_scale, column_colors, column_labels, genedata, samples, 
                         pvalues, title = "Rel. Expression", cluster_rows = TRUE, 
                         cluster_cols = TRUE, max_rows = 20, scale_by_row = TRUE) {
  
  annotation.row <- genedata %>%
    filter(ENSG %in% rownames(data)) 
  
  annotation.col <- samples %>%
    filter(sampleID %in% colnames(data)) %>%
    dplyr::select(c(class))
  
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
                  top_annotation = HeatmapAnnotation(class = anno_block(labels = column_labels,
                                                                        gp = gpar(fill = column_colors))),
                  heatmap_legend_param = list(legend_height = unit(10, "cm"),
                                              legend_width = unit(1.5, "cm"),
                                              title_position = "lefttop-rot"))
    
    draw(hm)
    
  } else {
    hm <- Heatmap(as.matrix(data), col = color_scale,  
                  name = title, cluster_rows = cluster_rows, column_split = annotation.col$class, 
                  show_row_dend = FALSE, show_column_names = FALSE, show_row_names = FALSE, 
                  row_title = NULL, column_title = NULL, cluster_columns = cluster_cols, 
                  top_annotation = HeatmapAnnotation(class = anno_block(labels = column_labels,
                                                                        gp = gpar(fill = column_colors))),
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
                  top_annotation = HeatmapAnnotation(class = anno_block(labels = column_labels,
                                                                        gp = gpar(fill = column_colors))),
                  heatmap_legend_param = list(legend_height = unit(10, "cm"),
                                              legend_width = unit(1.5, "cm"),
                                              title_position = "lefttop-rot"))
    
    draw(hm)
    
  }
}

#' Plot Volcano
#'
#' This function generates a volcano plot from DESeq2 results data, showing log2 fold change versus adjusted p-value.
#' It highlights significant genes and optionally caps log2 fold change values at a specified limit. It also adds gene labels
#' for the most significant genes.
#'
#' @param data DESeq2 results data.
#' @param title Title of the plot.
#' @param max_labels Maximum number of labels to show.
#' @param lfc_limit Limit for log fold change.
#' @return A ggplot object.
#' @export
plot_volcano <- function(data, title = NA, max_labels = NA, lfc_limit = NA) {
  data.plot <- data %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    mutate(genelabels = "", 
           out_of_bounds = (ifelse(is.na(vp_lfc_limit), 0, abs(log2FoldChange) > vp_lfc_limit) * sign(log2FoldChange)),
           log2FoldChange_capped = ifelse(out_of_bounds != 0, vp_lfc_limit * sign(log2FoldChange), log2FoldChange),
           padj = ifelse(padj == 0, .Machine$double.xmin, padj))
  
  if (is.na(max_labels) || sum(data.plot$threshold) < max_labels) {
    data.plot$genelabels[data.plot$threshold] <- data.plot$GeneSymbol[data.plot$threshold]
  } else if (max_labels > 0) {
    data.plot$genelabels[data.plot$threshold][1:max_labels] <- data.plot$GeneSymbol[data.plot$threshold][1:max_labels]
  }
  
  maxFC <- max(abs(data.plot$log2FoldChange))
  if (!is.na(lfc_limit) && lfc_limit < maxFC) {
    xlim <- lfc_limit
  } else {
    xlim <- maxFC * 1.04
  }
  
  breaks_y <- 10^-seq.int(0, -floor(log10(min(data.plot$padj))), by=20)
  
  
  ggplot(data.plot, aes(x = log2FoldChange_capped, y = padj)) +
    geom_point(data = subset(data.plot, out_of_bounds == 0), aes(colour=threshold), alpha = 0.5) +
    geom_point(data = subset(data.plot, out_of_bounds == -1), aes(colour=threshold), shape = "\u25c4", size=2) +
    geom_point(data = subset(data.plot, out_of_bounds == 1), aes(colour=threshold), shape = "\u25BA", size=2) +
    geom_hline(yintercept = 0.05, linetype = 2) +
    geom_vline(xintercept = c(-0.58, 0.58), linetype = 2) +
    geom_text_repel(aes(label = genelabels)) +
    scale_x_continuous(breaks = scales::pretty_breaks(), limits = c(-xlim, xlim), expand = expansion(0.01)) +
    scale_y_continuous(trans = c("log10", "reverse"), breaks = breaks_y) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "darkred")) +
    ggtitle(title) +
    xlab("log2 Fold Change") +
    ylab("adj. p-value") +
    theme_pubr() +
    theme(legend.position = "none")
}

