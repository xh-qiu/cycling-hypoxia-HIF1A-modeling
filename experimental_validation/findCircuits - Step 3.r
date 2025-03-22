#%%
# Load necessary libraries
library(data.table)
library(openxlsx)

# Define thresholds for p-value and log2 fold change
pval_threshold <- 0.05
l2fc_threshold <- 0.5

# Load DESeq2 results for WT and Normoxia
de_WT <- as.data.table(openxlsx::read.xlsx("../results/DESeq2_results_WT_with_pattern.xlsx", sheet = 1))
de_Normoxia <- as.data.table(openxlsx::read.xlsx("../results/DESeq2_results_Normoxia.xlsx", sheet = 1))

# Define an empty data table with columns TF, Target, Relationship
TFtarget_table <- data.table(TF = character(), Target = character(), Relationship = character())

# Define HIF1A to gene TF target relationships
up_HIF1A_genes <- de_WT[pvalue_Hypoxia_vs_Control_WT < 
pval_threshold & log2FoldChange_Hypoxia_vs_Control_WT > l2fc_threshold & !is.na(hgnc_symbol) & hgnc_symbol != "", .(Target = hgnc_symbol)]
down_HIF1A_genes <- de_WT[pvalue_Hypoxia_vs_Control_WT < 
pval_threshold & log2FoldChange_Hypoxia_vs_Control_WT < l2fc_threshold & !is.na(hgnc_symbol) & hgnc_symbol != "",.(Target = hgnc_symbol)]

#%%
# Add columns TF and Relationship to up_HIF1A_genes and down_HIF1A_genes
up_HIF1A_genes[, `:=` (TF = "HIF1A", Relationship = "up")]
down_HIF1A_genes[, `:=` (TF = "HIF1A", Relationship = "down")]

# Merge up and down HIF1A genes into the TFtarget_table
TFtarget_table <- rbind(TFtarget_table, up_HIF1A_genes, use.names = TRUE, fill = TRUE)
TFtarget_table <- rbind(TFtarget_table, down_HIF1A_genes, use.names = TRUE, fill = TRUE)

#%%

# Get p53 regulated genes
down_p53_genes <- de_Normoxia[pvalue.p53KO_vs_WT < 
pval_threshold & log2FoldChange.p53KO_vs_WT > l2fc_threshold & !is.na(hgnc_symbol) & hgnc_symbol != "", .(Target = hgnc_symbol)]
up_p53_genes <- de_Normoxia[pvalue.p53KO_vs_WT < 
pval_threshold & log2FoldChange.p53KO_vs_WT < l2fc_threshold & !is.na(hgnc_symbol) & hgnc_symbol != "", .(Target = hgnc_symbol)]

# Add columns TF and Relationship to up_p53_genes and down_p53_genes
up_p53_genes[, `:=` (TF = "TP53", Relationship = "up")]
down_p53_genes[, `:=` (TF = "TP53", Relationship = "down")]

# Merge up and down p53 genes into the TFtarget_table
TFtarget_table <- rbind(TFtarget_table, up_p53_genes, use.names = TRUE, fill = TRUE)
TFtarget_table <- rbind(TFtarget_table, down_p53_genes, use.names = TRUE, fill = TRUE)
#%%
# Get Notch1 regulated genes
down_notch1_genes <- de_Normoxia[pvalue.Notch1KO_vs_WT < 
pval_threshold & log2FoldChange.Notch1KO_vs_WT > l2fc_threshold & !is.na(hgnc_symbol) & hgnc_symbol != "", .(Target = hgnc_symbol)]
up_notch1_genes <- de_Normoxia[pvalue.Notch1KO_vs_WT < 
pval_threshold & log2FoldChange.Notch1KO_vs_WT < l2fc_threshold & !is.na(hgnc_symbol) & hgnc_symbol != "", .(Target = hgnc_symbol)]

# Add columns TF and Relationship to up_Notch1_genes and down_Notch1_genes
up_notch1_genes[, `:=` (TF = "NOTCH1", Relationship = "up")]
down_notch1_genes[, `:=` (TF = "NOTCH1", Relationship = "down")]

# Merge up and down Notch1 genes into the TFtarget_table
TFtarget_table <- rbind(TFtarget_table, up_notch1_genes, use.names = TRUE, fill = TRUE)
TFtarget_table <- rbind(TFtarget_table, down_notch1_genes, use.names = TRUE, fill = TRUE)
openxlsx::write.xlsx(TFtarget_table, "../results/TFtarget_table.xlsx", asTable = TRUE)

# Define a function to analyze specific FFLs
analyze_specific_FFLs <- function(data, X, Y) {
  results <- list()
  setDT(data)
  
  # Retrieve relationships
  X_targets <- data[TF == X, .(Target, Relationship)]
  Y_targets <- data[TF == Y, .(Target, Relationship)]
  X_Y_relation <- data[TF == X & Target == Y, .(Relationship)]
  sign_X_Y <- ifelse(nrow(X_Y_relation) > 0 && X_Y_relation$Relationship[1] == "up", "+", "-")
  
  # Check if X_Y_relation is empty and handle it
  if (nrow(X_Y_relation) == 0) {
    sign_X_Y <- ""  # Or handle as needed
  }

  # Merge targets
  common_targets <- merge(X_targets, Y_targets, by = "Target", suffixes = c(".X", ".Y"))

  for (i in seq_len(nrow(common_targets))) {
    row <- common_targets[i]
    Z <- row$Target
    sign_X_Z <- ifelse(row$Relationship.X == "up", "+", "-")
    sign_Y_Z <- ifelse(row$Relationship.Y == "up", "+", "-")
    signs <- paste(sign_X_Y, sign_Y_Z, sign_X_Z, sep="")
    
    # Assigning type based on the pattern of signs
    type <- switch(signs,
                   "+-+" = "Incoherent Type 1",
                   "---" = "Incoherent Type 2",
                   "++-" = "Incoherent Type 3",
                   "-++" = "Incoherent Type 4",
                   "+++" = "Coherent Type 1",
                   "-+-" = "Coherent Type 2",
                   "+--" = "Coherent Type 3",
                   "--+" = "Coherent Type 4",
                   "Other/Unmatched pattern")
    type_description <- switch(signs,
                   "+-+" = "X -> Y, Y -| Z, X -> Z",
                   "---" = "X -| Y, Y -| Z, X -| Z",
                   "++-" = "X -> Y, Y -> Z, X -| Z",
                   "-++" = "X -| Y, Y -> Z, X -> Z",
                   "+++" = "X -> Y, Y -> Z, X -> Z",
                   "-+-" = "X -| Y, Y -| Z, X -> Z",
                   "+--" = "X -> Y, Y -| Z, X -| Z",
                   "--+" = "X -| Y, Y -> Z, X -> Z",
                   "Other/Unmatched pattern")

    results[[length(results) + 1]] <- list(Type = type, X = X, Y = Y, Z = Z, Signs = signs, Description = type_description)
  }

  return(results)
}

# Convert TFtarget_table to data.table
setDT(TFtarget_table)

# Add relationships for TP53 and NOTCH1
TFtarget_table[(TF == "HIF1A") & (Target == "TP53"), Relationship := "up"]
TFtarget_table[(TF == "HIF1A") & (Target == "NOTCH1"), Relationship := "up"]

# Analyze FFLs
default_results_TP53 <- analyze_specific_FFLs(TFtarget_table, "HIF1A", "TP53")
default_results_NOTCH1 <- analyze_specific_FFLs(TFtarget_table, "HIF1A", "NOTCH1")

# Convert results to a data )table
results_dt_TP53 <- rbindlist(lapply(default_results_TP53, as.data.table), fill = TRUE)
results_dt_NOTCH1 <- rbindlist(lapply(default_results_NOTCH1, as.data.table), fill = TRUE)

results_dt <- rbind(results_dt_TP53, results_dt_NOTCH1, fill = TRUE)

# Save results to an Excel file
openxlsx::write.xlsx(results_dt, paste0("../results/IFFL_results.xlsx"), asTable = TRUE)

# Confirm the results have been saved
print("IFFL results have been successfully saved to 'IFFL_results.xlsx'")

