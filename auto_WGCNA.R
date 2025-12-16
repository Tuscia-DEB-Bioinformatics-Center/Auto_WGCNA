parse_arguments <- function() {
  option_list <- list(
    make_option(
      c("-d", "--workdir"),
      type = "character",
      help = "Working directory absolute path"
    ),
    make_option(
      c("-c", "--countmatrix"),
      type = "character",
      help = "Count matrix file absolute path"
    ),
    make_option(
      c("-p", "--phenodata"),
      type = "character",
      help = "Phenodata CSV file",
      default = NA
    ),
    make_option(
      c("-t", "--threads"),
      type = "integer",
      help = "Threads number",
      default = NA
    ),
    make_option(
      c("-m", "--mode"),
      type = "character",
      help = "Script execution mode [automatic (default) / manual]",
      default = "automatic"
    ),
    make_option(
      c("-y", "--threshold"),
      type = "double",
      help = "Threshold value",
      default = NA
    )
  )

  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)

  return(opt)
}

count_normalization <- function(dt_counts) {
  # Create DESeq2 object
  deseq_input <- as.matrix(dt_counts[, -1])
  row.names(deseq_input) <- dt_counts$gene_id

  col_data <- data.frame(row.names = colnames(deseq_input))

  dds <- DESeqDataSetFromMatrix(
    countData = deseq_input,
    colData = col_data,
    design = ~ 1
  )

  # Calculate normalization factors and apply VST
  dds <- estimateSizeFactors(dds)
  vsd <- vst(dds, blind = TRUE)      # blind true because no testing conditions
  vst_mat <- assay(vsd)

  # Traspose table for WGCNA: rows = samples, columns = genes
  vst_for_wgcna <- as.data.frame(t(vst_mat))

  # Export csv table
  write.table(
    vst_for_wgcna,
    "WGCNA_input_VST.csv",
    sep = ",",
    row.names = TRUE,
    quote = FALSE
  )

  zip("WGCNA_input_VST.zip", files = "WGCNA_input_VST.csv")

  file.remove("WGCNA_input_VST.csv")

  return(list(vst_for_wgcna = vst_for_wgcna, vst_mat = vst_mat))
}

threshold_calculation_manual <- function(vst_for_wgcna, max_power) {
  # Choose a set of soft-thresholding powers
  powers <- c(1, seq(2, max_power, by = 2))

  # Call the network topology analysis function
  sft <- pickSoftThreshold(vst_for_wgcna, powerVector = powers, verbose = 1)

  return(list(sft = sft, powers = powers))
}

threshold_calculation_automatic <- function(
  vst_for_wgcna,
  max_power,
  delta_min,
  min_length
) {
  # Choose a set of soft-thresholding powers
  powers <- c(1, seq(2, max_power, by = 2))

  # Call the network topology analysis function
  sft <- pickSoftThreshold(vst_for_wgcna, powerVector = powers, verbose = 1)
  fit <- sft$fitIndices[, 2]

  # Calculate absolute incremental differencies
  deltas <- abs(diff(fit))

  cat("delta --> ", deltas, "\n")

  plateaus <- deltas < delta_min

  # Identify TRUE consecutives sequences
  rle_plateau <- rle(plateaus)
  ends <- cumsum(rle_plateau$lengths)
  starts <- ends - rle_plateau$lengths + 1

  results <- list()

  for (i in seq_along(rle_plateau$values)) {
    if (rle_plateau$values[i] && rle_plateau$lengths[i] >= min_length) {
      indexes <- starts[i]:ends[i]

      results[[length(results) + 1]] <- data.frame(
        start_power = powers[indexes[1]],
        end_power = powers[ends[i] + 1],
        mean_value = mean(fit[indexes])
      )
    }
  }

  # Return max power if no plateau where found
  if (length(results) == 0) {
    cat("No plateau found, returning max power ...")
    return(list(sft = sft, picked_power = max(fit), powers = powers))
  }

  cat("Plateau found: \n")
  print(results)

  # Select longest plateau
  results_df <- do.call(rbind, results)
  lengths <- results_df$end_power - results_df$start_power
  longest <- which.max(lengths)
  picked_power <- results_df$mean_value[longest]

  cat("Picked power: \n")
  print(results_df$mean_value[longest])

  return(list(sft = sft, picked_power = picked_power, powers = powers))
}

plot_threshold_charts <- function(sft, picked_power, powers) {
  # Save plot as TIFF image
  tiff(
    "charts/threshold.tiff",
    width = 5760,
    height = 3240,
    units = "px",
    res = 600,
    compression = "zip"
  )

  par(mfrow = c(1, 2))
  cex1 <- 0.8

  plot(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    xlab = "Soft Threshold (power)",
    ylab = "Scale Free Topology Model Fit, signed R^2",
    main = paste("Scale independence")
  )
  text(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    labels = powers,
    cex = cex1,
    col = "red"
  )
  if (picked_power > 0) {
    abline(h = picked_power, col = "red")
  }
  plot(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    xlab = "Soft Threshold (power)",
    ylab = "Mean Connectivity",
    type = "n",
    main = paste("Mean connectivity")
  )
  text(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    labels = powers,
    cex = cex1, col = "red"
  )

  dev.off()

  # Save plot as PNG image
  png(
    "charts/threshold.png",
    width = 960,
    height = 540,
    units = "px",
    res = 150
  )

  par(mfrow = c(1, 2))
  cex1 <- 0.8

  plot(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    xlab = "Soft Threshold (power)",
    ylab = "Scale Free Topology Model Fit, signed R^2",
    main = paste("Scale independence")
  )
  text(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    labels = powers,
    cex = cex1,
    col = "red"
  )
  if (picked_power > 0) {
    abline(h = picked_power, col = "red")
  }
  plot(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    xlab = "Soft Threshold (power)",
    ylab = "Mean Connectivity",
    type = "n",
    main = paste("Mean connectivity")
  )
  text(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    labels = powers,
    cex = cex1, col = "red"
  )

  dev.off()
}

network_build <- function(vst_for_wgcna, picked_power) {
  temp_cor <- cor
  cor <- WGCNA::cor         # Force it to use WGCNA cor function 

  network <- blockwiseModules(
    vst_for_wgcna,                       # count matrix
    power = picked_power * 10,           # soft threshold
    networkType = "signed",              # network type: signed = more restrictive correlation, unsigned = less restrictive correlation         
    deepSplit = 2,                       # 0 - 4, 0 = (less modules with bigger size), 4 = (more modules with smaller size)
    pamRespectsDendro = FALSE,           # consent structure changes
    # detectCutHeight = 0.75,
    minModuleSize = 30,                  # min number of genes to build a module
    maxBlockSize = 4000,                 # max size of genes block
    reassignThreshold = 0,               # genes reassignment between modules during analisys (default = 1e-6, 0 = no reassignment)
    mergeCutHeight = 0.25,               # threshold to merge similar modules (default = 0.25)
    saveTOMs = FALSE,                     # archive the run results in TOM file (saves time)
    numericLabels = TRUE,                # modules label format (TRUE = numbers, FALSE = colors)
    verbose = 1
  )

  return(network)
}

plot_dendogram <- function(network) {
  # Convert labels to colors for plotting
  colors <- labels2colors(network$colors)

  # Save plot as TIFF image
  tiff(
    "charts/dendogram.tiff",
    width = 5760,
    height = 3240,
    units = "px",
    res = 600,
    compression = "zip"
  )
  # Plot the dendrogram and the module colors underneath
  plotDendroAndColors(
    network$dendrograms[[1]],             # dendogram object to be visualized
    colors[network$blockGenes[[1]]],      # colors vector
    "Module colors",
    dendroLabels = FALSE,                 # NULL = display branch labels, FALSE = no display branch labels
    hang = 0.03,                          # branch's length (default = 0.1)
    addGuide = TRUE,                      # display color's guide
    guideHang = 0.05
  )
  dev.off()

  # Save plot as PNG image
  png(
    "charts/dendogram.png",
    width = 960,
    height = 540,
    units = "px",
    res = 150
  )
  # Plot the dendrogram and the module colors underneath
  plotDendroAndColors(
    network$dendrograms[[1]],             # dendogram object to be visualized
    colors[network$blockGenes[[1]]],      # colors vector
    "Module colors",
    dendroLabels = FALSE,                 # NULL = display branch labels, FALSE = no display branch labels
    hang = 0.03,                          # branch's length (default = 0.1)
    addGuide = TRUE,                      # display color's guide
    guideHang = 0.05
  )
  dev.off()

  return(colors)
}

save_colors_list <- function(network) {
  # Get the list of modules
  gene_modules <- data.frame(
    gid = names(network$colors),
    color = labels2colors(network$colors)
  )

  # Save list of colors to file
  cat("Saving modules list to file ... \n")
  write.table(
    gene_modules,
    file = "gene_modules.tsv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  cat("Modules list saved to file\n")

  return(gene_modules)
}

identify_eigengenes <- function(vst_for_wgcna, colors) {
  # Identify eigengenes gene of each module
  eigengenes_matrix <- moduleEigengenes(vst_for_wgcna, colors)$eigengenes

  # Order modules
  eigengenes_matrix <- orderMEs(eigengenes_matrix)
  module_order <- gsub("ME", "", names(eigengenes_matrix))

  # ...
  eigengenes_matrix$sample <- row.names(eigengenes_matrix)

  # Convert dataframe in sample, name, value
  eigengenes_matrix_m <- mutate(
    pivot_longer(eigengenes_matrix, -sample),
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

  # Save dataframe to TSV file
  write.table(
    eigengenes_matrix_m,
    file = "eigengenes_matrix_m.tsv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  return(
    list(
      eigengenes_matrix = eigengenes_matrix,
      eigengenes_matrix_m = eigengenes_matrix_m
    )
  )
}

calculate_modpheno_correlation <- function(pheno_dt, eigengenes_matrix) {
  # Convert non numeric data in numeric
  pheno_num_dt <- data.frame(
    lapply(pheno_dt, function(x) as.numeric(as.factor(x)))
  )

  rownames(pheno_num_dt) <- rownames(pheno_dt)

  # Create module x phenotype matrix
  module_phenodata_cor <- cor(
    eigengenes_matrix[, !(names(eigengenes_matrix) %in% c("sample"))],
    pheno_num_dt,
    use = "p"
  )

  module_pheno_pvalue <- corPvalueStudent(
    module_phenodata_cor,
    nSamples = nrow(eigengenes_matrix)
  )

  mpc_save <- cbind(module = rownames(module_phenodata_cor), module_phenodata_cor)
  rownames(mpc_save) <- NULL

  # Save dataframe to TSV file
  write.table(
    mpc_save,
    file = "module_phenodata_cor.tsv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  return(module_phenodata_cor)
}

plot_modpheno_heatmap <- function(
  module_phenodata_cor,
  pheno_dt,
  eigengenes_matrix
) {
  print(length(names(pheno_dt)))
  print(length(names(eigengenes_matrix)))

  # Save plot as TIFF image
  tiff(
    "charts/module_pheno_heatmap.tiff",
    width = 2 * length(names(pheno_dt)),
    height = 0.2 * length(names(eigengenes_matrix)),
    units = "in",
    res = 600,
    compression = "zip"
  )
  labeledHeatmap(
    Matrix = module_phenodata_cor,
    xLabels = names(pheno_dt),
    yLabels = names(eigengenes_matrix[, !(names(eigengenes_matrix) %in% c("sample"))]),
    ySymbols = names(eigengenes_matrix[, !(names(eigengenes_matrix) %in% c("sample"))]),
    colorLabels = FALSE,
    colors = blueWhiteRed(50),
    textMatrix = signif(module_phenodata_cor, 2),
    setStdMargins = FALSE,
    cex.lab = 0.3,
    cex.text = 0.2,
    zlim = c(-1, 1),
    main = "Module-Phenodata Correlations"
  )
  dev.off()

  # Save plot as PNG image
  png(
    "charts/module_pheno_heatmap.png",
    width = 2 * length(names(pheno_dt)),
    height = 0.2 * length(names(eigengenes_matrix)),
    units = "in",
    res = 150
  )
  labeledHeatmap(
    Matrix = module_phenodata_cor,
    xLabels = names(pheno_dt),
    yLabels = names(eigengenes_matrix[, !(names(eigengenes_matrix) %in% c("sample"))]),
    ySymbols = names(eigengenes_matrix[, !(names(eigengenes_matrix) %in% c("sample"))]),
    colorLabels = FALSE,
    colors = blueWhiteRed(50),
    textMatrix = signif(module_phenodata_cor, 2),
    setStdMargins = FALSE,
    cex.lab = 0.3,
    cex.text = 0.2,
    zlim = c(-1, 1),
    main = "Module-Phenodata Correlations"
  )
  dev.off()
}

plot_heatmap <- function(eigengenes_matrix_m, eigengenes_matrix) {
  heatmap <- ggplot(
    eigengenes_matrix_m,
    aes(x = sample, y = name, fill = value)
  ) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1, 1)) +
    theme(axis.text.x = element_text(angle = 90)
  ) +
  labs(title = "Module-trait Relationships", y = "Modules", fill = "correlation")

  # Export heatmap to TIFF file
  ggsave(
    "charts/heatmap.tiff",
    plot = heatmap,
    width = 0.5 * length(unique(eigengenes_matrix_m$sample)),
    height = 0.5 * length(unique(eigengenes_matrix_m$name)),
    dpi = 600,
    limitsize = FALSE,
    device = "tiff",
    compression = "zip"
  )

  # Export heatmap to PNG file
  ggsave(
    "charts/heatmap.png",
    plot = heatmap,
    width = 0.5 * length(unique(eigengenes_matrix_m$sample)),
    height = 0.5 * length(unique(eigengenes_matrix_m$name)),
    dpi = 150,
    limitsize = FALSE
  )
}

plot_expression_profiles <- function(vst_mat, gene_modules, workdir) {
  # Get the list of submodules (colors)
  colors <- unique(gene_modules$color)

  for (col in colors) {
    mask <- gene_modules$color %in% col

    submod <- subset(gene_modules, mask)

    row.names(gene_modules) <- gene_modules$gid

    subexpr <- vst_mat[submod$gid, ]

    submod_df <- data.frame(subexpr)
    submod_df$gid <- row.names(submod_df)

    # Make submod data frame in long format
    submod_df_long <- pivot_longer(submod_df, -gid)
    submod_df_long$module <- gene_modules[submod_df_long$gid, ]$color
    submod_df <- submod_df_long

    # Create modules profiles chart
    profiles <- ggplot(
      submod_df,
      aes(x = name, y = value, group = gid)
    ) + 
    geom_line(aes(color = module), alpha = 0.2) +
    scale_color_identity() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    facet_grid(rows = vars(module)) +
    labs(x = "condition", y = "normalized expression")

    # Create profile plot name
    profile_name_tiff <- paste("profile_", col, ".tiff", sep = "")
    profile_name_png <- paste("profile_", col, ".png", sep = "")

    # Export profiles to TIFF file
    ggsave(
      profile_name_tiff,
      plot = profiles,
      path = "module_profiles",
      width = 0.3 * length(unique(submod_df$name)),
      height = 4,
      dpi = 600,
      device = "tiff",
      compression = "zip"
    )

    # Export profiles to PNG file
    ggsave(
      profile_name_png,
      plot = profiles,
      path = "module_profiles",
      width = 0.3 * length(unique(submod_df$name)),
      height = 4,
      dpi = 150
    )
  }

  zip(
    "charts/module_profiles.zip",
    files = "module_profiles/"
  )

  unlink("module_profiles", recursive = TRUE)
}

plot_boxplots <- function(pheno_dt, colors, eigengenes_matrix_m) {
  pheno_df <- data.frame(sample = rownames(pheno_dt))
  pheno_df <- cbind(pheno_df, pheno_dt)

  table_df <- left_join(
    eigengenes_matrix_m,
    pheno_df,
    by = "sample"
  )

  write.table(
    table_df,
    file = "table.tsv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  for (col in unique(colors)) {
    mask <- table_df$name %in% col

    submod <- subset(table_df, mask)

    for (con in colnames(table_df)[4:ncol(table_df)]) {
      # Create boxplot chart
      boxplot <- ggplot(
        submod,
        aes(x = submod[[con]], y = value, fill = submod[[con]])
      ) +
        geom_boxplot(outlier.shape = NA, alpha = 0.6) +
        geom_jitter(width = 0.2, size = 1, alpha = 0.8) +
        geom_text_repel(
          aes(label = sample),
          size = 1,
          max.overlaps = Inf,
          position = position_jitter(width = 0.2),
          show.legend = FALSE
        ) +
        facet_wrap(~ name, scales = "free_y") +
        theme_bw() +
        labs(
          x = "Condition",
          y = "Normalized expression (Module Eigengene)",
          title = "Module expression for condition"
        ) +
        theme(
          legend.position = "none",
          strip.text = element_text(face = "bold")
        )

      # Create box plot name
      boxplot_tiff <- paste("boxplot_", col, "_", con, ".tiff", sep = "")
      boxplot_png <- paste("profile_", col, "_", con, ".png", sep = "")

      # Export boxplot to TIFF file
      ggsave(
        boxplot_tiff,
        plot = boxplot,
        path = "boxplots",
        width = 3 * length(unique(submod[[con]])),
        height = 4,
        dpi = 600,
        device = "tiff",
        compression = "zip"
      )

      # Export profiles to PNG file
      ggsave(
        boxplot_png,
        plot = boxplot,
        path = "boxplots",
        width = 3 * length(unique(submod[[con]])),
        height = 4,
        dpi = 150
      )
    }
  }

  zip(
    "charts/boxplots.zip",
    files = "boxplots/"
  )

  unlink("boxplots", recursive = TRUE)
}

identify_hub_genes <- function(vst_for_wgcna, colors, picked_power) {
  # Calculate hub genes for each module
  hubs <- chooseTopHubInEachModule(
    datExpr = vst_for_wgcna,
    colorh = colors,          # modules colors
    power = picked_power,     # power used for adjacency network
    type = "signed"           # entered network type
  )

  # Print result to stdout
  print(hubs)

  # Convert result in data frame
  hubs_df <- data.frame(Module = names(hubs), TopHub = hubs)

  # Save result to files
  write.table(hubs_df,
    file = "hub_genes.tsv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

plot_top_gene_network <- function(
  vst_for_wgcna,
  picked_power,
  network,
  colors
) {
  # Identify top genes and order in descend order
  sel_genes <- names(
    sort(
      apply(vst_for_wgcna, 2, var),
      decreasing = TRUE
    )[1:1000]
  )

  # Calculate adjacency matrix
  adjacency_matrix <- adjacency(
    vst_for_wgcna[, sel_genes],
    power = picked_power,
    type = "signed"
  )

  # Calculate Topological Overlap Matrix
  tom <- TOMsimilarity(adjacency_matrix)
  colnames(tom) <- sel_genes
  rownames(tom) <- sel_genes

  # Keep top 20% of connections
  tom_trimmed <- quantile(tom[upper.tri(tom)], 0.80)
  cat("TOM threshold (80%):", tom_trimmed, "\n")

  adjacency_threshold <- tom
  adjacency_threshold[adjacency_threshold < tom_trimmed] <- 0

  # Create graph object
  network_graph <- graph_from_adjacency_matrix(
    adjacency_threshold,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

  # Remove unconnected nodes
  network_graph <- delete_vertices(
    network_graph,
    degree(network_graph) == 0
  )

  # Get module colors for these genes
  node_colors <- colors[match(sel_genes, names(network$colors))]

  # Calculate node sizes based on connectivity
  node_sizes <- log10(degree(network_graph) + 1) * 3 + 3

  # Save plot as TIFF image
  tiff(
    "charts/coexpression_network.tiff",
    width = 10,
    height = 7,
    units = "in",
    res = 600,
    compression = "zip"
  )
  par(mar = c(0, 0, 3, 0))
  set.seed(123)
  plot(
    network_graph,
    vertex.size = node_sizes,
    vertex.label = NA,
    vertex.color = node_colors,
    vertex.frame.color = NA,
    edge.width = 0.2,
    edge.color = "gray8",
    layout = layout_with_fr(network_graph),
    main = "Network\nTop 1000 genes"
  )

  unique_colors <- unique(node_colors[node_colors != "grey"])

  legend(
    "topright",
    legend = unique_colors,
    col = unique_colors,
    pch = 19,
    cex = 1,
    title = "Modules",
    bty = "n"
  )
  dev.off()
}

write_exec_time <- function(
  start_time,
  end_time,
  deseq_time_start,
  deseq_time_end,
  threshold_calc_time_start,
  threshold_calc_time_end,
  net_build_time_start,
  net_build_time_end,
  exp_plot_time_start,
  exp_plot_time_end
) {

  total_exec_time <- difftime(end_time, start_time, units = "secs")
  deseq_time <- difftime(deseq_time_end, deseq_time_start, units = "secs")
  threshold_time <- difftime(threshold_calc_time_end, threshold_calc_time_start, units = "secs")
  net_build_time <- difftime(net_build_time_end, net_build_time_start, units = "secs")
  exp_plot_time <- difftime(exp_plot_time_end, exp_plot_time_start, units = "secs")

  cat(
    "Total execution time:", total_exec_time, "seconds", "\n",
    "Deseq2 execution time:", deseq_time, "seconds", "\n",
    "Threshold calculation time:", threshold_time, "seconds", "\n",
    "Net build execution time:", net_build_time, "seconds", "\n",
    "Expression profiles plot time:", exp_plot_time, "seconds", "\n",
    "\n", "-------------------------------------------------", "\n\n",
    file = "execution_time.txt",
    append = TRUE
  )
}

main <- function() {
  # Load libraries
  suppressWarnings(suppressPackageStartupMessages({
    library(WGCNA)
    library(tidyverse)
    library(ggrepel)
    library(DESeq2)
    library(optparse)
    library(igraph)
  }))

  # Setting string not as factor
  options(stringAsFactors = FALSE)

  # Parse command arguments
  opt <- parse_arguments()

  if (!is.na(opt$threads)) {
    if (opt$threads > 1) {
      # Enable multithreading with specified threads
      enableWGCNAThreads(nThreads = opt$threads)
    } else {
      stop("Threads number must be an integer > 1")
    }
  } else {
    # Enable multithreading with max-1 threads
    enableWGCNAThreads()
  }

  # Set current working directory
  setwd(opt$workdir)
  cat("Working directory settled \n")

  # Only if in automatic mode or first part of manual mode
  if (is.na(opt$threshold)) {
    # Create directory
    dir.create("charts", showWarnings = TRUE, recursive = FALSE, mode = "0777")
    dir.create("module_profiles", showWarnings = TRUE, recursive = FALSE, mode = "0777")
    dir.create("boxplots", showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }

  # Reading the raw data (rows are the sample and columns the genes)
  dt_counts <- read.csv(opt$countmatrix, sep = ",")
  cat("count matrix read \n")

  # Reading the phenodata data (columns are sample, phenodata1, phenodata2, ...)
  if (!is.na(opt$phenodata)) {
    pheno_dt <- read.csv(opt$phenodata, sep = "\t")

    # Add first column a row names
    row.names(pheno_dt) <- pheno_dt[, 1]

    # Remove first column of phenodata dataframe
    pheno_dt <- pheno_dt[, -1]
  }

  start_time <- Sys.time()

  # Count normalization (using DESeq2)
  deseq_time_start <- Sys.time()

  count_normalization_res <- count_normalization(dt_counts)
  vst_for_wgcna <- count_normalization_res$vst_for_wgcna
  vst_mat <- count_normalization_res$vst_mat

  deseq_time_end <- Sys.time()
  cat("Count normalization done \n")

  # Threshold calculation
  cat("threshold calculation started \n")
  threshold_calc_time_start <- Sys.time()

  # Soft threshold calculation
  if (opt$mode == "automatic") {
    sft_results <- threshold_calculation_automatic(
      vst_for_wgcna = vst_for_wgcna,
      max_power = 50,
      delta_min = 0.01,
      min_length = 3
    )

    sft <- sft_results$sft
    picked_power <- sft_results$picked_power
    powers <- sft_results$powers
  } else {
    if (opt$mode == "manual" && is.na(opt$threshold)) {
      sft_results <- threshold_calculation_manual(
        vst_for_wgcna = vst_for_wgcna,
        max_power = 50
      )

      sft <- sft_results$sft
      powers <- sft_results$powers

      # Plot Softhreshold charts
      plot_threshold_charts(sft, 0, powers)

      # Save data
      save(sft, powers, file = "sft_p.RData")

      cat(
        "Check the 'threshold' chart and choose a value for soft threshold.", "\n",
        "Restart the script selection with option -m manual and selected threshold value in option -y \n"
      )
      return()
    }

    if (opt$mode == "manual" && !is.na(opt$threshold)) {
      # Load previous saved sft and power vector data
      load("sft_p.RData")

      # Use user's selected power
      picked_power <- opt$threshold
    }
  }

  threshold_calc_time_end <- Sys.time()
  cat("threshold calculation finished \n")

  # Plot Softhreshold charts
  cat("plot threshold started \n")
  plot_threshold_charts(sft, picked_power, powers)
  cat("plot threshold finished \n")

  # Network build
  cat("network build started \n")
  net_build_time_start <- Sys.time()

  network <- network_build(vst_for_wgcna, picked_power)

  net_build_time_end <- Sys.time()
  cat("network build finished \n")

  # Cluster plot
  cat("plot dendogram started \n")
  colors <- plot_dendogram(network)
  cat("plot dendogram finished \n")

  # Save gene modules
  gene_modules <- save_colors_list(network)

  # Identify and save eigengenes for each module
  eigen_results <- identify_eigengenes(vst_for_wgcna, colors)
  eigengenes_matrix_m <- eigen_results$eigengenes_matrix_m
  eigengenes_matrix <- eigen_results$eigengenes_matrix

  # Calculate and plot module-phenodata correlations
  if (!is.na(opt$phenodata)) {
    # Calculate module phenotype correlation
    module_phenodata_cor <- calculate_modpheno_correlation(
      pheno_dt,
      eigengenes_matrix
    )

    # Plot and save module x phenotype correlation heatmap chart
    plot_modpheno_heatmap(
      module_phenodata_cor,
      pheno_dt,
      eigengenes_matrix
    )
  }

  # Generate and save heatmap chart
  plot_heatmap(eigengenes_matrix_m, eigengenes_matrix)

  # Generate and save gene modules expression profiles
  cat("Generating gene expression profiles charts ... \n")
  exp_plot_time_start <- Sys.time()

  plot_expression_profiles(vst_mat, gene_modules, opt$workdir)

  exp_plot_time_end <- Sys.time()

  cat("Generating box plots ...")
  plot_boxplots(pheno_dt, colors, eigengenes_matrix_m)

  # Identify hub gene for each module
  cat("Identifying hub gene for each module ... \n")
  identify_hub_genes(vst_for_wgcna, colors, picked_power)
  cat("Hub genes correctly identified.")

  # Plot top gene network graph
  plot_top_gene_network(vst_for_wgcna, picked_power, network, colors)

  end_time <- Sys.time()

  # Write execution time to file
  write_exec_time(
    start_time,
    end_time,
    deseq_time_start,
    deseq_time_end,
    threshold_calc_time_start,
    threshold_calc_time_end,
    net_build_time_start,
    net_build_time_end,
    exp_plot_time_start,
    exp_plot_time_end
  )

}

main()