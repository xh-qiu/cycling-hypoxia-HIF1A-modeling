#%%
library(data.table)
library(ggplot2)
library(biomaRt)

#%%
wt_de <- openxlsx::read.xlsx("../results/DESeq2_results_WT.xlsx", sheet = 1)
wt_de <- data.table(wt_de)

p53_de <- openxlsx::read.xlsx("../results/DESeq2_results_p53KO.xlsx", sheet = 1)
p53_de <- data.table(p53_de)

Notch1_de <- openxlsx::read.xlsx("../results/DESeq2_results_Notch1KO.xlsx", sheet = 1)
Notch1_de <- data.table(Notch1_de)

#%%n define the pattern function
# Function to determine the pattern based on fold changes and p-values
# Function to determine the ordering pattern for samples H, N, and O.
getPattern <- function(l2FC_H_vs_N, l2FC_O_vs_N, l2FC_O_vs_H,
                       pval_H_vs_N, pval_O_vs_N, pval_O_vs_H,
                       pval_threshold = 0.5, l2fc_threshold = 1) {
  
  # convert NA to 0
  l2FC_H_vs_N <- ifelse(is.na(l2FC_H_vs_N), 0, l2FC_H_vs_N)
  l2FC_O_vs_N <- ifelse(is.na(l2FC_O_vs_N), 0, l2FC_O_vs_N)
  l2FC_O_vs_H <- ifelse(is.na(l2FC_O_vs_H), 0, l2FC_O_vs_H)

  # convert pvals NA to 1
  pval_H_vs_N <- ifelse(is.na(pval_H_vs_N), 1, pval_H_vs_N)
  pval_O_vs_N <- ifelse(is.na(pval_O_vs_N), 1, pval_O_vs_N)
  pval_O_vs_H <- ifelse(is.na(pval_O_vs_H), 1, pval_O_vs_H)
  
  # Helper function to compare a pair of samples.
  # It returns:
  #   "eq" if the p-value is above threshold OR the absolute log2 fold-change is too small,
  #   "gt" if the log2 fold-change is significantly positive (i.e. first sample > second sample),
  #   "lt" if the log2 fold-change is significantly negative (i.e. first sample < second sample).
  comparePair <- function(l2fc, pval, pval_threshold, l2fc_threshold) {
    if (pval > pval_threshold || abs(l2fc) < l2fc_threshold) {
      return("eq")
    } else if (l2fc > 0) {
      return("gt")
    } else {
      return("lt")
    }
  }
  
  # Determine the relation between each pair.
  # For H vs N: a positive l2FC_H_vs_N means H > N.
  rel_HN <- comparePair(l2FC_H_vs_N, pval_H_vs_N, pval_threshold, l2fc_threshold)
  # For O vs N: a positive l2FC_O_vs_N means O > N.
  rel_ON <- comparePair(l2FC_O_vs_N, pval_O_vs_N, pval_threshold, l2fc_threshold)
  # For O vs H: a positive l2FC_O_vs_H means O > H.
  rel_OH <- comparePair(l2FC_O_vs_H, pval_O_vs_H, pval_threshold, l2fc_threshold)
  
  # (Optional) If all three pairwise comparisons are approximate, then all three samples are similar.
  if (rel_HN == "eq" && rel_ON == "eq" && rel_OH == "eq") {
    # Not one of the allowed outputs but we flag it.
    return("H~N~O")
  }
  
  # ---------------------------
  # Handle cases with one pair approximately equal
  # ---------------------------
  
  # Case 1: H and N are approximately equal.
  if (rel_HN == "eq") {
    # Check how O compares to N (or H; they are ~equal).
    if (rel_ON == "eq") {
      # If O and N are also ~equal, then all three are approximately equal.
      return("H~N~O")
    } else if (rel_ON == "gt") {
      # O > N, so the (H, N) group is lower than O.
      # Use the option "N~H<O" (order within the equal pair is irrelevant).
      return("N~H<O")
    } else if (rel_ON == "lt") {
      # O < N, so ordering is O < (H, N). We choose "O<N~H".
      return("O<N~H")
    }
  }
  
  # Case 2: O and N are approximately equal.
  if (rel_ON == "eq") {
    # Compare H with N (or O; they are ~equal).
    if (rel_HN == "gt") {
      # H > N so (O, N) is lower than H.
      # Use "O~N<H" (order within the equal pair is irrelevant).
      return("O~N<H")
    } else if (rel_HN == "lt") {
      # H < N so ordering is H < (O, N), i.e. "H<O~N".
      return("H<O~N")
    }
  }
  
  # Case 3: O and H are approximately equal.
  if (rel_OH == "eq") {
    # Compare N with H (or O, since they are ~equal).
    # We use the H vs N comparison (rel_HN).
    if (rel_HN == "gt") {
      # H > N so the group (O, H) is higher than N.
      # Use "N<O~H" to indicate N is lowest.
      return("N<O~H")
    } else if (rel_HN == "lt") {
      # H < N so the group (O, H) is lower than N.
      # Use "O~H<N" to indicate N is highest.
      return("O~H<N")
    }
  }
  
  # ---------------------------
  # Handle strict comparisons (no approximate equalities)
  # ---------------------------
  
  # At this point, none of the pairwise comparisons are approximate ("eq").
  # We have strict relations for each pair. The six allowed orderings (always written
  # in ascending order, i.e. from lowest to highest) are:
  #   "O<N<H", "N<O<H", "N<H<O", "O<H<N", "H<O<N", "H<N<O"
  
  # Use the three strict relations to deduce the overall order.
  if (rel_HN == "gt" && rel_ON == "gt" && rel_OH == "gt") {
    # H > N, O > N, and O > H, so the order is: N < H < O.
    return("N<H<O")
  }
  if (rel_HN == "gt" && rel_ON == "gt" && rel_OH == "lt") {
    # H > N, O > N, but O < H implies: N < O < H.
    return("N<O<H")
  }
  if (rel_HN == "lt" && rel_ON == "lt" && rel_OH == "gt") {
    # H < N, O < N, but O > H implies: H < O < N.
    return("H<O<N")
  }
  if (rel_HN == "lt" && rel_ON == "lt" && rel_OH == "lt") {
    # H < N, O < N, and O < H implies: O < H < N.
    return("O<H<N")
  }
  if (rel_HN == "gt" && rel_ON == "lt") {
    # H > N and O < N gives: O < N < H.
    return("O<N<H")
  }
  if (rel_HN == "lt" && rel_ON == "gt") {
    # H < N and O > N gives: H < N < O.
    return("H<N<O")
  }
  
  # If none of the above conditions are met, the result is ambiguous.
  return("Impossible")
}



#%%
oscPattern <- function(l2FC_H_vs_N, l2FC_O_vs_N, l2FC_O_vs_H,
                       pval_H_vs_N, pval_O_vs_N, pval_O_vs_H,
                       pval_threshold = 0.05, l2fc_threshold = 0.2, l2fc_threshold_between = 0.05) {

  # convert NA pvals to 1
  pval_H_vs_N <- ifelse(is.na(pval_H_vs_N), 1, pval_H_vs_N)
  pval_O_vs_N <- ifelse(is.na(pval_O_vs_N), 1, pval_O_vs_N)
  pval_O_vs_H <- ifelse(is.na(pval_O_vs_H), 1, pval_O_vs_H)

  # convert NA l2FCs to 0
  l2FC_H_vs_N <- ifelse(is.na(l2FC_H_vs_N), 0, l2FC_H_vs_N)
  l2FC_O_vs_N <- ifelse(is.na(l2FC_O_vs_N), 0, l2FC_O_vs_N)
  l2FC_O_vs_H <- ifelse(is.na(l2FC_O_vs_H), 0, l2FC_O_vs_H)

  if ( (pval_O_vs_N < pval_threshold) & (pval_O_vs_H < pval_threshold)) { # O is differentially expressed vs both N and H
    if ( (l2FC_O_vs_N > l2fc_threshold) & (l2FC_O_vs_H > l2fc_threshold) ) { # O is upregulated vs both N and H
      if (l2FC_H_vs_N > 0) { # H is upregulated vs N
        return("OscUpExtreme")
      } else {
        return("OscUpOpposite")
      }
    } else if ( (l2FC_O_vs_N < -l2fc_threshold) & (l2FC_O_vs_H < -l2fc_threshold) ) { # O is downregulated vs both N and H
      if (l2FC_H_vs_N < 0) { # H is downregulated vs N
        return("OscDnExtreme")
      } else {
        return("OscDnOpposite")
      }
    }
  }

  # check if O is between H and N
  if ( (l2FC_O_vs_N > l2fc_threshold_between) & (l2FC_O_vs_H < -l2fc_threshold_between) ) { # O is upregulated vs N and downregulated vs H
    return("Between")
  }
  if ( (l2FC_O_vs_N < -l2fc_threshold_between) & (l2FC_O_vs_H > l2fc_threshold_between) ) { # O is downregulated vs N and upregulated vs H
    return("Between")
  }

  return("Other")
}



#%% do wt

# Apply the updated function to all rows with p-values
wt_de[, pattern := mapply(getPattern, 
                          log2FoldChange_Hypoxia_vs_Control_WT, 
                          log2FoldChange_Oscillation_vs_Control_WT, 
                          log2FoldChange_Oscillation_vs_Hypoxia_WT,
                          pvalue_Hypoxia_vs_Control_WT,
                          pvalue_Oscillation_vs_Control_WT,
                          pvalue_Oscillation_vs_Hypoxia_WT,
                          l2fc_threshold = 0.25, pval_threshold = 0.01)]

wt_de[, type_osc := mapply(oscPattern, 
                              log2FoldChange_Hypoxia_vs_Control_WT, 
                              log2FoldChange_Oscillation_vs_Control_WT, 
                              log2FoldChange_Oscillation_vs_Hypoxia_WT,
                              pvalue_Hypoxia_vs_Control_WT,
                              pvalue_Oscillation_vs_Control_WT,
                              pvalue_Oscillation_vs_Hypoxia_WT)]


# Save the updated table to a new Excel file
openxlsx::write.xlsx(wt_de, "../results/DESeq2_results_WT_with_pattern.xlsx", asTable = TRUE)

# Check the patterns
table(wt_de$pattern)
table(wt_de$type_osc)
#%% do p53ko
p53_de[, pattern := mapply(getPattern, 
                           log2FoldChange_Hypoxia_vs_Control_p53KO, 
                           log2FoldChange_Oscillation_vs_Control_p53KO, 
                           log2FoldChange_Oscillation_vs_Hypoxia_p53KO,
                           pvalue_Hypoxia_vs_Control_p53KO,
                           pvalue_Oscillation_vs_Control_p53KO,
                           pvalue_Oscillation_vs_Hypoxia_p53KO,
                           l2fc_threshold = 0.05, pval_threshold = 0.05)]

p53_de[, type_osc := mapply(oscPattern, 
                              log2FoldChange_Hypoxia_vs_Control_p53KO, 
                              log2FoldChange_Oscillation_vs_Control_p53KO, 
                              log2FoldChange_Oscillation_vs_Hypoxia_p53KO,
                              pvalue_Hypoxia_vs_Control_p53KO,
                              pvalue_Oscillation_vs_Control_p53KO,
                              pvalue_Oscillation_vs_Hypoxia_p53KO)]


openxlsx::write.xlsx(p53_de, "../results/DESeq2_results_p53KO_with_pattern.xlsx", asTable = TRUE)

table(p53_de$pattern)

#%% do Notch1ko
Notch1_de[, pattern := mapply(getPattern, 
                              log2FoldChange_Hypoxia_vs_Control_Notch1KO, 
                              log2FoldChange_Oscillation_vs_Control_Notch1KO, 
                              log2FoldChange_Oscillation_vs_Hypoxia_Notch1KO,
                              pvalue_Hypoxia_vs_Control_Notch1KO,
                              pvalue_Oscillation_vs_Control_Notch1KO,
                              pvalue_Oscillation_vs_Hypoxia_Notch1KO,
                              l2fc_threshold = 0.05, pval_threshold = 0.05)]

Notch1_de[, type_osc := mapply(oscPattern, 
                              log2FoldChange_Hypoxia_vs_Control_Notch1KO, 
                              log2FoldChange_Oscillation_vs_Control_Notch1KO, 
                              log2FoldChange_Oscillation_vs_Hypoxia_Notch1KO,
                              pvalue_Hypoxia_vs_Control_Notch1KO,
                              pvalue_Oscillation_vs_Control_Notch1KO,
                              pvalue_Oscillation_vs_Hypoxia_Notch1KO)]

openxlsx::write.xlsx(Notch1_de, "../results/DESeq2_results_Notch1KO_with_pattern.xlsx", asTable = TRUE)

table(Notch1_de$pattern)

#%%
# merge the tables, putting a suffix on column name pattern to keep track of the source, and also remove all columns that contain the word "Lactate"

wt_de <- openxlsx::read.xlsx("../results/DESeq2_results_WT_with_pattern.xlsx", sheet = 1)
p53_de <- openxlsx::read.xlsx("../results/DESeq2_results_p53KO_with_pattern.xlsx", sheet = 1)
Notch1_de <- openxlsx::read.xlsx("../results/DESeq2_results_Notch1KO_with_pattern.xlsx", sheet = 1)

# convert all to data tables
wt_de <- data.table(wt_de)
p53_de <- data.table(p53_de)
Notch1_de <- data.table(Notch1_de)

# Remove columns with "Lactate"
wt_de <- wt_de[, !grepl("Lactate", names(wt_de)), with = FALSE]
p53_de <- p53_de[, !grepl("Lactate", names(p53_de)), with = FALSE]
Notch1_de <- Notch1_de[, !grepl("Lactate", names(Notch1_de)), with = FALSE]

# Add a suffix to the pattern column name to keep track of the source
setnames(wt_de, "pattern", "pattern_WT")
setnames(p53_de, "pattern", "pattern_p53KO")
setnames(Notch1_de, "pattern", "pattern_Notch1KO")

setnames(wt_de, "type_osc", "type_osc_WT")
setnames(p53_de, "type_osc", "type_osc_p53KO")
setnames(Notch1_de, "type_osc", "type_osc_Notch1KO")

# Merge the tables
merged_de <- merge(wt_de, p53_de, by = c("ensembl_gene_id", "hgnc_symbol"))
merged_de <- merge(merged_de, Notch1_de, by = c("ensembl_gene_id", "hgnc_symbol"))

# Save the merged table to a new Excel file
openxlsx::write.xlsx(merged_de, "../results/DESeq2_results_merged_with_pattern.xlsx", asTable = TRUE, overwrite = TRUE)

# add up the occurrences of each pattern, and make a table with columns pattern, WT, p53KO, Notch1KO. 
# The numbers in the table should be the counts of each pattern in each dataset.

pattern_counts <- data.table(table(merged_de$pattern_WT, merged_de$pattern_p53KO, merged_de$pattern_Notch1KO))
setnames(pattern_counts, c("V1", "V2", "V3"), c("pattern_WT", "pattern_p53KO", "pattern_Notch1KO"))

pattern_counts_osc <- data.table(table(merged_de$type_osc_WT, merged_de$type_osc_p53KO, merged_de$type_osc_Notch1KO))
setnames(pattern_counts_osc, c("V1", "V2", "V3"), c("type_osc_WT", "type_osc_p53KO", "type_osc_Notch1KO"))

# Save the pattern counts table to a new Excel file
openxlsx::write.xlsx(list(pattern_counts = pattern_counts, type_counts = pattern_counts_osc), "../results/pattern_counts.xlsx", asTable = TRUE)
    

    

#%%
# download GO gene sets from MSigDB
# https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

# read gmt gene sets file "../../../downloaded/c5.go.bp.v2024.1.Hs.symbols.gmt"
read_gmt_file <- function(file_path) {
  # Read the file
  lines <- readLines(file_path)
  
  # Initialize an empty list to store the pathways
  pathways <- list()
  
  # Process each line
  for (line in lines) {
    # Split the line by tabs
    elements <- unlist(strsplit(line, "\t"))
    
    # The first element is the pathway name
    pathway_name <- elements[1]
    
    # The second element is the pathway description
    pathway_description <- elements[2]
    
    # The remaining elements are the genes
    genes <- elements[3:length(elements)]
    
    # Add the genes to the pathways list with the pathway name as the key
    pathways[[pathway_name]] <- genes
    
    # Store the description as an attribute
    attr(pathways[[pathway_name]], "description") <- pathway_description
  }
  
  return(pathways)
}

# Usage:
go_pathways <- read_gmt_file("../data/downloaded/c5.go.bp.v2024.1.Hs.symbols.gmt")

kegg_pathways <- read_gmt_file("../data/downloaded/c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt")

#%%
library(fgsea)

# get the list of genes with patterns "HNO" in the WT dataset
wt_de <- openxlsx::read.xlsx("../results/DESeq2_results_WT_with_pattern.xlsx", sheet = 1)
wt_de <- data.table(wt_de)

# Get the genes with the specified patterns
genes_HNO <- wt_de[hgnc_symbol %in% wt_de$hgnc_symbol[wt_de$pattern == "H<N<O"], hgnc_symbol]
genes_HNO <- genes_HNO[!duplicated(genes_HNO)]
genes_HNO <- genes_HNO[!is.na(genes_HNO)]
genes_HNO <- genes_HNO[!genes_HNO == ""]

genes_NHO <- wt_de[hgnc_symbol %in% wt_de$hgnc_symbol[wt_de$pattern == "N<H<O"], hgnc_symbol]
genes_NHO <- genes_NHO[!duplicated(genes_NHO)]
genes_NHO <- genes_NHO[!is.na(genes_NHO)]
genes_NHO <- genes_NHO[!genes_NHO == ""]


genes_OHN <- wt_de[hgnc_symbol %in% wt_de$hgnc_symbol[wt_de$pattern == "O<H<N"], hgnc_symbol]
genes_OHN <- genes_OHN[!duplicated(genes_OHN)]
genes_OHN <- genes_OHN[!is.na(genes_OHN)]
genes_OHN <- genes_OHN[!genes_OHN == ""]


genes_ONH <- wt_de[hgnc_symbol %in% wt_de$hgnc_symbol[wt_de$pattern == "O<N<H"], hgnc_symbol]
genes_ONH <- genes_ONH[!duplicated(genes_ONH)]
genes_ONH <- genes_ONH[!is.na(genes_ONH)]
genes_ONH <- genes_ONH[!genes_ONH == ""]

univ <- unique(wt_de$hgnc_symbol)
univ <- univ[!duplicated(univ)]
univ <- univ[!is.na(univ)]
univ <- univ[!univ == ""]

go_HNO = fgsea::fora(go_pathways, genes_HNO, universe = univ)
go_NHO = fgsea::fora(go_pathways, genes_NHO, universe = univ)
go_OHN = fgsea::fora(go_pathways, genes_OHN, universe = univ)
go_ONH = fgsea::fora(go_pathways, genes_ONH, universe = univ)

kegg_HNO = fgsea::fora(kegg_pathways, genes_HNO, universe = univ)
kegg_NHO = fgsea::fora(kegg_pathways, genes_NHO, universe = univ)
kegg_OHN = fgsea::fora(kegg_pathways, genes_OHN, universe = univ)
kegg_ONH = fgsea::fora(kegg_pathways, genes_ONH, universe = univ)

# save all the GO results to an Excel file
openxlsx::write.xlsx(list(go_HNO = go_HNO, go_NHO = go_NHO, go_OHN = go_OHN, go_ONH = go_ONH), "../results/go_oscillation_wt.xlsx", asTable = TRUE)

# save all the KEGG results to an Excel file
openxlsx::write.xlsx(list(kegg_HNO = kegg_HNO, kegg_NHO = kegg_NHO, kegg_OHN = kegg_OHN, kegg_ONH = kegg_ONH), "../results/kegg_oscillation_wt.xlsx", asTable = TRUE)



#%%
# plan to do:
## Plan 01 : find genes that are upregulated in WT with hypoxia (hypoxia vs control), but not up-regulated in p53KO with hypoxia (hypoxia vs control)
library(data.table)
library(openxlsx)

# Load data
wt_de <- openxlsx::read.xlsx("../results/DESeq2_results_WT_with_pattern.xlsx", sheet = 1)
wt_de <- data.table(wt_de)

p53_de <- openxlsx::read.xlsx("../results/DESeq2_results_p53KO_with_pattern.xlsx", sheet = 1)
p53_de <- data.table(p53_de)

Notch1_de <- openxlsx::read.xlsx("../results/DESeq2_results_Notch1KO_with_pattern.xlsx", sheet = 1)
Notch1_de <- data.table(Notch1_de)

# Define thresholds for determining upregulation
logFC_threshold <- 0  
padj_threshold <- 0.05

# Filter WT genes that are upregulated under hypoxia vs control
wt_upregulated <- wt_de[
  log2FoldChange_Hypoxia_vs_Control_WT > logFC_threshold & 
  padj_Hypoxia_vs_Control_WT < padj_threshold]

# Filter p53KO genes that are upregulated under hypoxia vs control
p53ko_upregulated <- p53_de[
  log2FoldChange_Hypoxia_vs_Control_p53KO > logFC_threshold & 
  padj_Hypoxia_vs_Control_p53KO < padj_threshold]


# Filter Notch1KO genes that are upregulated under hypoxia vs control 
Notch1_upregulated <- Notch1_de[
  log2FoldChange_Hypoxia_vs_Control_Notch1KO > logFC_threshold & 
  padj_Hypoxia_vs_Control_Notch1KO < padj_threshold]


# Identify genes upregulated in WT but not in p53KO
specific_up_wt_not_p53ko <- wt_upregulated[!ensembl_gene_id %in% p53ko_upregulated$ensembl_gene_id]

# Output these genes
openxlsx::write.xlsx(specific_up_wt_not_p53ko, "../results/Specific_Upregulated_in_WT_not_p53KO.xlsx", asTable = TRUE)

 
#%%


## Plan 02 :  find genes that are upregulated in WT with hypoxia (hypoxia vs control), but not up-regulated in Notch1KO with hypoxia (hypoxia vs control)
specific_up_wt_not_Notch1ko <- wt_upregulated[!ensembl_gene_id %in% Notch1_upregulated$ensembl_gene_id]
# Output these genes
openxlsx::write.xlsx(specific_up_wt_not_Notch1ko, "../results/Specific_Upregulated_in_WT_not_Notch1KO.xlsx", asTable = TRUE)

## Plan 03 :  find genes that are downregulated in WT with hypoxia (hypoxia vs control), but not down-regulated in p53KO with hypoxia (hypoxia vs control)
specific_down_wt_not_p53ko <- wt_de[
  log2FoldChange_Hypoxia_vs_Control_WT < -logFC_threshold & 
  padj_Hypoxia_vs_Control_WT < padj_threshold][
    !ensembl_gene_id %in% p53_de[
      log2FoldChange_Hypoxia_vs_Control_p53KO < -logFC_threshold & 
      padj_Hypoxia_vs_Control_p53KO < padj_threshold]$ensembl_gene_id]
# Output these genes
openxlsx::write.xlsx(specific_down_wt_not_p53ko, "../results/Specific_Downregulated_in_WT_not_p53KO.xlsx", asTable = TRUE)


## Plan 04 :  find genes that are downregulated in WT with hypoxia (hypoxia vs control), but not down-regulated in Notch1KO with hypoxia (hypoxia vs control)
specific_down_wt_not_Notch1ko <- wt_de[
  log2FoldChange_Hypoxia_vs_Control_WT < -logFC_threshold & 
  padj_Hypoxia_vs_Control_WT < padj_threshold][
    !ensembl_gene_id %in% Notch1_de[
      log2FoldChange_Hypoxia_vs_Control_Notch1KO < -logFC_threshold & 
      padj_Hypoxia_vs_Control_Notch1KO < padj_threshold]$ensembl_gene_id]
# Output these genes
openxlsx::write.xlsx(specific_down_wt_not_Notch1ko, "../results/Specific_Downregulated_in_WT_not_Notch1KO.xlsx", asTable = TRUE)

#%%
## Plan 05 :  do pathway analysis on these genes, separately, and save it in Excel files
# Load the gene sets
go_pathways <- read_gmt_file("../data/downloaded/c5.go.bp.v2024.1.Hs.symbols.gmt") # GO gene sets
kegg_pathways <- read_gmt_file("../data/downloaded/c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt") # KEGG gene sets  
# analysis on these genes in Specific_Upregulated_in_WT_not_p53KO.xlsx  
# analysis on these genes in Specific_Upregulated_in_WT_not_Notch1KO.xlsx

#%%

# function to remove empty and NA values from a vector
remove_empty_na <- function(x) {
  x <- x[!is.na(x)]
  x <- x[x != ""]
  return(x)
}

up_not_p53 <- fgsea::fora(go_pathways, remove_empty_na(specific_up_wt_not_p53ko$hgnc_symbol), universe = univ)
dn_not_p53 <- fgsea::fora(go_pathways, remove_empty_na(specific_down_wt_not_p53ko$hgnc_symbol), universe = univ)



up_not_Notch1 <- fgsea::fora(go_pathways, remove_empty_na(specific_up_wt_not_Notch1ko$hgnc_symbol), universe = univ)
dn_not_Notch1 <- fgsea::fora(go_pathways, remove_empty_na(specific_down_wt_not_Notch1ko$hgnc_symbol), universe = univ)

# save all the GO results to an Excel file
openxlsx::write.xlsx(list(up_not_p53 = up_not_p53, dn_not_p53 = dn_not_p53, up_not_Notch1 = up_not_Notch1, dn_not_Notch1 = dn_not_Notch1), "../results/go_hypoxia_KO.xlsx", asTable = TRUE)


#%%
# Plan 06 :  Venn diagrams showing the number of oscillation patterns in WT, p53KO, and Notch1KO
# Load necessary libraries
library(VennDiagram)
library(grid)
library(data.table)

# Convert merged_de to a data.table
setDT(merged_de)

# Function to Extract Gene Sets by Pattern
get_gene_sets_by_pattern <- function(pattern) {
  list(
    WT = merged_de[pattern_WT == pattern, ensembl_gene_id],
    p53KO = merged_de[pattern_p53KO == pattern, ensembl_gene_id],
    Notch1KO = merged_de[pattern_Notch1KO == pattern, ensembl_gene_id]
  )
}

# Define All Patterns to Analyze
patterns <- c("N<H<O", "O<N<H", "O<H<N",  "H<N<O")

# Loop Through Patterns and Generate Venn Diagrams
for (pattern in patterns) {
  # Get gene sets for this pattern
  venn_list <- get_gene_sets_by_pattern(pattern)

  # Output file name
  output_file <- paste0("../figures/venn_", pattern, ".png")
  
  # Generate and Save Venn Diagram for Each Pattern
  png(output_file, width = 2000, height = 2000, res = 300)
  
  venn_plot <- venn.diagram(
    x = venn_list, 
    filename = NULL,  # Prevents auto-saving
    category.names = c("WT", "p53KO", "Notch1KO"),  # Labels
    fill = c("cornflowerblue", "green", "yellow"),  # Colors
    alpha = 0.50, 
    cex = 1.5, 
    cat.cex = 1.5,
    cat.pos = c(-45, 45, 180),  # Adjust label positioning
    cat.dist = c(0.05, 0.07, 0.05),  # Adjust spacing
    margin = 0.15,  # Margin around the Venn diagram
    label.col = "black",  # Count color
    print.mode = "raw",  # Show counts directly
    fontface = "bold"  # Make numbers bold
  )

  # Draw the Venn Diagram
  grid.newpage()
  grid.draw(venn_plot)

  # ✅ Only add pattern label at the top (remove Notch1KO duplication)
  #grid.text(pattern, x = 0.50, y = 0.92, gp = gpar(fontsize = 20, col = "black", fontface = "bold"))
  grid.text(paste0("Oscillation Patterns Venn Diagram - ", pattern), 
          x = 0.50, y = 0.92, 
          gp = gpar(fontsize = 20, col = "black", fontface = "bold"))

  # Close the PNG device
  dev.off()

  # Confirm file creation
  if (file.exists(output_file)) {
    print(paste("✅ Venn diagram for", pattern, "saved as:", output_file))
  } else {
    stop(paste("❌ Failed to save the Venn diagram for", pattern))
  }
}

#%% Create a combined bar chart for all patterns
library(ggplot2)
library(data.table)

# Initialize an empty list to store counts for all patterns
all_pattern_counts <- list()

patterns <- c("N<H<O", "O<N<H", "O<H<N", "H<N<O")

# Collect data for all patterns
for (pattern in patterns) {
  genes_in_patterns <- get_gene_sets_by_pattern(pattern)
  wt_genes <- genes_in_patterns$WT
  p53_genes <- intersect(wt_genes, genes_in_patterns$p53KO) 
  Notch1_genes <- intersect(wt_genes, genes_in_patterns$Notch1KO)
  
  all_pattern_counts[[pattern]] <- data.table(
    group = c("WT", "p53KO", "Notch1KO"),
    count = c(length(wt_genes), length(p53_genes), length(Notch1_genes)),
    pattern = pattern
  )
}
all_pattern_counts
# Combine all data
combined_counts <- rbindlist(all_pattern_counts)

combined_counts$group <- factor(combined_counts$group, levels = c("WT", "p53KO", "Notch1KO"))

# Create the combined bar chart
combined_chart <- ggplot(combined_counts, aes(x = pattern, y = count, fill = group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  scale_fill_manual(values = c("WT" = "cornflowerblue", 
                              "p53KO" = "forestgreen", 
                              "Notch1KO" = "goldenrod")) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.text = element_text(color = "black")
  ) +
  labs(
    title = "Gene Counts by Pattern",
    x = "Pattern",
    y = "Number of Genes",
    fill = "Group"
  )
print(combined_chart)
# Save the combined bar chart
ggsave("../figures/oscillation_patterns_combined_bar_chart.tiff", 
       combined_chart, 
       width = 5, 
       height = 3, 
       dpi = 300,
       compression = "lzw")







