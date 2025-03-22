library(data.table)
library(openxlsx)
library(fgsea)
library(ggplot2)

# Load WT differential expression results
wt_de <- data.table(read.xlsx("../results/DESeq2_results_WT_with_pattern.xlsx"))

# Load GO pathways
go_pathways <- gmtPathways("../data/downloaded/c5.go.bp.v2024.1.Hs.symbols.gmt")

# Patterns to analyze
patterns <- c("N<H<O", "O<H<N", "O<N<H", "H<N<O")
# Initialize a list to store results
go_top5_results_list <- list()  # Stores top 5 pathways for each pattern
go_all_results_list <- list()   # Stores all pathways for each pattern

all_genes = wt_de$hgnc_symbol
all_genes = all_genes[!is.na(all_genes) & all_genes != ""]
# Run fgsea separately for each pattern and create plots
for (this_pattern in patterns) {

    # Select genes belonging to the current pattern
    genes_pattern <- wt_de[(pattern == this_pattern) & !is.na(hgnc_symbol) & hgnc_symbol != "", hgnc_symbol]
    

    fora_res = fgsea::fora(genes = genes_pattern, pathways = go_pathways, universe = all_genes, minSize = 10, maxSize = 200)
    fora_res[, enrichment := overlap/size]
    # Select top 5 pathways by smallest p-value
    top5_go <- fora_res[order(pval)][1:5]
    
    # Select top 5 pathways by smallest p-value
    # Store results in the list
    go_top5_results_list[[this_pattern]] <- top5_go  # Store top 5 pathways
    go_all_results_list[[this_pattern]] <- fora_res[order(pval)]  # Sort all pathways by p-value

    # Create barplot
    plot <- ggplot(top5_go, aes(x = reorder(pathway, enrichment), y = enrichment, fill = enrichment)) +
        geom_bar(stat = "identity", width = 0.6, color = "black")+
        scale_fill_distiller(palette = "YlGnBu", direction = 1) +
        coord_flip() +
        labs(title = paste("Top 5 GO pathways for pattern", this_pattern),
             x = "GO Pathways",
             y = "Enrichment (overlap/size)",
             fill = "Enrichment (overlap/size)") +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1),
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              axis.text = element_text(size = 10, color = "black"),
              legend.title = element_text(size = 12),
              legend.text = element_text(size = 10),
              plot.title = element_text(size = 14, hjust = 0.5))

    # Save plot as PNG file
    output_filename <- paste0("../figures/Top5_GO_Pathways_WT_", gsub("[<>]", "", this_pattern), ".tiff")
    ggsave(output_filename, plot = plot, width = 4+max(nchar(top5_go$pathway))/10, height = 2, 
        device = "tiff", compression = "lzw+p", units = "in")
    
    cat("Plot successfully saved as", output_filename, "\n")
}

# Save Top 5 GO pathways results in an Excel file with separate sheets
write.xlsx(go_top5_results_list, "../results/Top5_GO_Pathways_WT.xlsx", asTable = TRUE, overwrite = TRUE)

# Save All GO pathways results in an Excel file with separate sheets
write.xlsx(go_all_results_list, "../results/All_GO_Pathways_WT.xlsx", asTable = TRUE, overwrite = TRUE)

cat("All results saved successfully.\n")

