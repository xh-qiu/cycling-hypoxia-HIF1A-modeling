# Load required libraries
library(data.table)
library(openxlsx)
library(ggplot2)
library(reshape2)

# Read TPM data
tpm_values <- as.data.table(openxlsx::read.xlsx("../results/tpm_mat.xlsx"))

# Read pattern assignment file
pattern_data <- as.data.table(openxlsx::read.xlsx("../results/DESeq2_results_merged_with_pattern.xlsx"))

# Define WT sample columns
wt_columns <- c("WT_Control_1", "WT_Control_2", "WT_Control_3",
                "WT_Hypoxia_1", "WT_Hypoxia_2", "WT_Hypoxia_3",
                "WT_Oscillation_1", "WT_Oscillation_2", "WT_Oscillation_3")

# Merge TPM values with pattern assignment using 'ensembl_gene_id'
merged_data <- merge(tpm_values, pattern_data[, .(ensembl_gene_id, pattern_WT)], by = "ensembl_gene_id")

# Define function to create heatmaps for each pattern using ggplot2
makeHeatMap <- function(pattern, tpm_data, dozscore = TRUE) {
    
    # Filter data for the selected pattern
    tpm_subset <- tpm_data[pattern_WT == pattern, ..wt_columns]
    
    # Check if any genes match the pattern
    if (nrow(tpm_subset) == 0) {
        warning(paste("No matching genes found for pattern:", pattern))
        return(NULL)
    }
    
    # Add gene names explicitly
    tpm_subset$Gene <- tpm_data[pattern_WT == pattern, ensembl_gene_id]

    # remove rows with all zero values
    tpm_subset <- tpm_subset[rowSums(tpm_subset[, -"Gene", with = FALSE]) != 0]

    # Reorder columns to place Gene first
    setcolorder(tpm_subset, c("Gene", wt_columns))

    # Perform Z-score normalization per row independently if requested
    if (dozscore) {
        tpm_subset[, (wt_columns) := {
            # Calculate row means and standard deviations
            row_mean = rowMeans(.SD)
            row_sd = apply(.SD, 1, sd)
            # Apply z-score normalization for each column
            lapply(.SD, function(x) (x - row_mean) / row_sd)
        }, .SDcols = wt_columns]
    }
    
    # Convert to long format for ggplot2
    heatmap_data <- melt(tpm_subset, id.vars = "Gene", variable.name = "Sample", value.name = "z-score")
    # Create heatmap using ggplot2
    p <- ggplot(heatmap_data, aes(x = Sample, y = Gene, fill = `z-score`)) +
        geom_tile() +
        scale_fill_gradient2(low = "#0017b0", mid = "white", high = "#741e1e", midpoint = 0, limits = c(-3, 3)) + 
        theme_classic() +
        labs(title = paste("Heatmap for", pattern), x = "Sample", y = "Gene", fill = "Z-score") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.text.y = element_blank())
    
    # Save heatmap to file
    file_path <- paste0("../figures/heatmaps/heatmap_", pattern, ".png")
    ggsave(file_path, plot = p, width = 6, height = 2 + nrow(tpm_subset)/300, units = "in", dpi = 300)
}


# List of patterns to generate heatmaps for
patterns <- unique(pattern_data$pattern_WT)

# Generate heatmaps for each pattern
for (pattern in patterns) {
    cat("Processing pattern:", pattern, "\n")
    makeHeatMap(pattern, merged_data, dozscore = TRUE)
}
