#%% load libraries
library(data.table)
library(ggplot2)

results_patterns <- as.data.table(openxlsx::read.xlsx("../results/DESeq2_results_merged_with_pattern.xlsx"))
IFFL_results <- as.data.table(openxlsx::read.xlsx("../results/IFFL_results.xlsx"))


head(results_patterns$pattern_WT)

# function to remove < symbol from strings
remove_less_than <- function(x) {
    gsub("<", "", x)
}

#%% find Osc genes with pattern

Osc_types <- c("OscUpExtreme", "OscUpOpposite", "OscDnExtreme", "OscDnOpposite")
Non_Osc_types <- c("Between")
wt_osc_genes <- results_patterns[type_osc_WT %in% Osc_types, hgnc_symbol]
p53_osc_genes <- results_patterns[type_osc_p53KO %in% Osc_types, hgnc_symbol]
notch1_osc_genes <- results_patterns[type_osc_Notch1KO %in% Osc_types, hgnc_symbol]

wt_non_osc_genes <- results_patterns[type_osc_WT %in% Non_Osc_types, hgnc_symbol]
p53_non_osc_genes <- results_patterns[type_osc_p53KO %in% Non_Osc_types, hgnc_symbol]
notch1_non_osc_genes <- results_patterns[type_osc_Notch1KO %in% Non_Osc_types, hgnc_symbol]

#%% 

# Create categorization for oscillation status
IFFL_results[, WT_Osc := fcase(
    Z %in% wt_osc_genes , "Osc",
    Z %in% wt_non_osc_genes , "Non-Osc",
    default = "Other"
)]

IFFL_results[, p53_Osc := fcase(
    Z %in% p53_osc_genes , "Osc",
    Z %in% p53_non_osc_genes , "Non-Osc",
    default = "Other"
)]

IFFL_results[, Notch1_Osc := fcase(
    Z %in% notch1_osc_genes , "Osc",
    Z %in% notch1_non_osc_genes , "Non-Osc",
    default = "Other"
)]

# Merge patterns with IFFL_results instead of doing row-by-row updates
IFFL_results <- merge(IFFL_results, 
                     results_patterns[, .(hgnc_symbol, type_osc_WT, type_osc_p53KO, type_osc_Notch1KO)],
                     by.x = "Z", 
                     by.y = "hgnc_symbol",
                     all.x = TRUE)

#%%
IFFL_results[grepl(Type, pattern = "Incoherent")][WT_Osc == "Osc"][Y == "TP53"][p53_Osc == "Non-Osc"] 

table(IFFL_results$Y)

# save this in excel
openxlsx::write.xlsx(IFFL_results[grepl(Type, pattern = "Incoherent")][WT_Osc == "Osc"][Y == "TP53"][p53_Osc == "Non-Osc"], 
                    "../results/WT_OSC_TP53.xlsx", asTable = TRUE)


IFFL_results[grepl(Type, pattern = "Incoherent")][WT_Osc == "Osc"][Y == "NOTCH1"][Notch1_Osc == "Non-Osc"]

# save this in excel
openxlsx::write.xlsx(IFFL_results[grepl(Type, pattern = "Incoherent")][WT_Osc == "Osc"][Y == "NOTCH1"][Notch1_Osc == "Non-Osc"], 
                    "../results/WT_OSC_NOTCH1.xlsx", asTable = TRUE)

#%% read in tpm values
tpm_values <- as.data.table(openxlsx::read.xlsx("../results/tpm_mat.xlsx"))


makeHeatMap <- function(hgnc_list = NULL, tpm_ = tpm_values, dozscore = TRUE, dodendo = TRUE,
                        genotypes = c("WT", "p53KO", "Notch1KO"),
                        conditions = c("Control", "Hypoxia", "Lactate", "Oscillation")) {
    
    # Filter genes and handle empty input
    hgnc_list = hgnc_list[!is.na(hgnc_list)]
    if (length(hgnc_list) == 0) {
        stop("No genes provided in hgnc_list")
    }
    
    tpm_subset = tpm_[hgnc_symbol %in% hgnc_list]
    if (nrow(tpm_subset) == 0) {
        stop("No matching genes found in TPM data")
    }
    
    # Create pattern to match desired columns
    col_pattern = paste0(
        "^(", 
        paste(genotypes, collapse="|"),
        ")_(",
        paste(conditions, collapse="|"),
        ")_[0-9]+$"
    )
    
    # Select columns matching pattern plus hgnc_symbol
    cols_to_keep = c("hgnc_symbol", grep(col_pattern, colnames(tpm_subset), value = TRUE))
    dt_this = tpm_subset[, ..cols_to_keep]
    
    # Convert to matrix for processing
    mat_this <- as.matrix(dt_this[, -1])
    rownames(mat_this) <- dt_this$hgnc_symbol
    mat_this <- mat_this[genefilter::rowSds(mat_this) > 1e-6, , drop = FALSE]
    
    if (dozscore) {
        # Handle single row case for mean
        if (nrow(mat_this) == 1) {
            mat_this_colnames <- colnames(mat_this)
            mat_this_rownames <- rownames(mat_this)
            mat_this = matrix(mat_this - mean(mat_this), nrow=1)
            colnames(mat_this) <- mat_this_colnames
            rownames(mat_this) <- mat_this_rownames
        } else {
            mat_this = sweep(x = mat_this, MARGIN = 1, FUN = "-", STATS = rowMeans(mat_this))
        }
        
        # Handle single row case for sd
        sds <- if (nrow(mat_this) == 1) {
            sd(as.numeric(mat_this))
        } else {
            genefilter::rowSds(mat_this)
        }
        
        # Add check for zero standard deviation
        sds[sds < 1e-10] <- 1
        
        # Apply scaling
        if (nrow(mat_this) == 1) {
            mat_this_colnames <- colnames(mat_this)
            mat_this_rownames <- rownames(mat_this)
            mat_this = matrix(as.numeric(mat_this) / sds, nrow=1)
            colnames(mat_this) <- mat_this_colnames
            rownames(mat_this) <- mat_this_rownames
        } else {
            mat_this = sweep(x = mat_this, MARGIN = 1, FUN = "/", STATS = sds)
        }
    }
    
    # Prepare data for plotting
    z_df = reshape2::melt(mat_this)
    if (dozscore) {
        names(z_df) <- c("Gene", "Sample", "z-score")
    } else {
        names(z_df) <- c("Gene", "Sample", "TPM")
    }
    
    # Add dendrogram ordering if requested
    if (dodendo && nrow(mat_this) > 1) {
        z.dendro <- as.dendrogram(hclust(d = dist(x = mat_this)))
        plot_dend = ggdendro::ggdendrogram(data = z.dendro, rotate = TRUE) + 
            theme(axis.text.y = element_text(size = 6), axis.text.x = element_blank())
        z_df$Gene <- factor(x = z_df$Gene, 
                          levels = rownames(mat_this)[order.dendrogram(z.dendro)], 
                          ordered = TRUE)
    } else {
        # For single gene or when dendrogram is not requested
        z_df$Gene = factor(x = as.character(z_df$Gene), levels = rev(hgnc_list))
    }
    
    # Extract condition and genotype from sample names
    z_df$genotype = sub("^([^_]+)_.*$", "\\1", z_df$Sample)
    z_df$condition = sub("^[^_]+_([^_]+)_.*$", "\\1", z_df$Sample)

    # order the genotypes and conditions
    z_df$genotype <- factor(z_df$genotype, levels = genotypes)
    z_df$condition <- factor(z_df$condition, levels = conditions)

    # find absolute max z-score
    abs_max_z <- max(abs(z_df$`z-score`))

    # Create the plot
    p <- ggplot(data = z_df, 
           mapping = aes(y = Gene, x = Sample, fill = `z-score`)) + 
        geom_tile() + 
        scale_fill_distiller(palette = "BrBG", limits = c(-abs_max_z, abs_max_z)) + 
        theme_minimal() +
        facet_grid(~ genotype + condition, scales = "free_x", space = "free_x") + 
        scale_x_discrete(position = "bottom") + 
        labs(y = "") + 
        theme(legend.position = "top", 
              axis.text = element_text(color = "black"), 
              axis.text.x = element_blank(),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.grid = element_blank(),
              strip.text = element_text(face = "bold"))
    
    return(p)
}

p53_circuit_osc <- IFFL_results[grepl(Type, pattern = "Incoherent")][WT_Osc == "Osc"][Y == "TP53"][p53_Osc == "Non-Osc"]

num_genes <- length(p53_circuit_osc$Z)
tiff(paste0("../figures/heatmaps/p53_circuit_together.tiff"), 
            width = 6, height = 1.5 + num_genes * 0.2, units = "in", res = 300, compression = "lzw+p")
makeHeatMap(hgnc_list = p53_circuit_osc$Z, tpm_ = tpm_values, dozscore = TRUE, dodendo = TRUE,
                genotypes = c("WT", "p53KO"),
                conditions = c("Control", "Hypoxia", "Oscillation"))
dev.off()


for (WT_pat in unique(p53_circuit_osc$type_osc_WT)) {
    p53_circuit_osc_genes <- IFFL_results[grepl(Type, pattern = "Incoherent")][WT_Osc == "Osc"][Y == "TP53"][p53_Osc == "Non-Osc"][type_osc_WT == WT_pat, Z]
    num_genes <- length(p53_circuit_osc_genes)
    # write num_genes and pattern to stdout
    cat(paste0("Number of genes: ", num_genes, "\nPattern: ", WT_pat, "\n"))
    tiff(paste0("../figures/heatmaps/p53_circuit_osc_genes_", remove_less_than(WT_pat), ".tiff"), 
            width = 6, height = 1.5 + num_genes * 0.2, units = "in", res = 300, compression = "lzw+p")
    print(makeHeatMap(hgnc_list = p53_circuit_osc_genes, tpm_ = tpm_values, dozscore = TRUE, dodendo = TRUE,
                genotypes = c("WT", "p53KO"),
                conditions = c("Control", "Hypoxia", "Oscillation")))
    dev.off()
}

notch1_circuit_osc <- IFFL_results[grepl(Type, pattern = "Incoherent")][WT_Osc == "Osc"][Y == "NOTCH1"][Notch1_Osc == "Non-Osc"]
num_genes <- length(notch1_circuit_osc$Z)
tiff(paste0("../figures/heatmaps/notch1_circuit_together.tiff"), 
            width = 6, height = 1.5 + num_genes * 0.2, units = "in", res = 300, compression = "lzw+p")
makeHeatMap(hgnc_list = notch1_circuit_osc$Z, tpm_ = tpm_values, dozscore = TRUE, dodendo = TRUE,
                genotypes = c("WT", "Notch1KO"),
                conditions = c("Control", "Hypoxia", "Oscillation"))
dev.off()

for (WT_pat in unique(notch1_circuit_osc$type_osc_WT)) {
    notch1_circuit_osc_genes <- IFFL_results[grepl(Type, pattern = "Incoherent")][WT_Osc == "Osc"][Y == "NOTCH1"][Notch1_Osc == "Non-Osc"][type_osc_WT == WT_pat, Z]
    num_genes <- length(notch1_circuit_osc_genes)
    cat(paste0("Number of genes: ", num_genes, "\nPattern: ", WT_pat, "\n"))
    tiff(paste0("../figures/heatmaps/notch1_circuit_osc_genes_", remove_less_than(WT_pat), ".tiff"), 
            width = 6, height = 1.5 + num_genes * 0.2, units = "in", res = 300, compression = "lzw+p")
    print(makeHeatMap(hgnc_list = notch1_circuit_osc_genes, tpm_ = tpm_values, dozscore = TRUE, dodendo = TRUE,
                genotypes = c("WT", "Notch1KO"),
                conditions = c("Control", "Hypoxia", "Oscillation")))
    dev.off()
}


