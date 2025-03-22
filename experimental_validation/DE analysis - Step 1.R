#%%
library(ggplot2)
library(DESeq2)
library(data.table)
library(biomaRt)
library(ggrepel)

#%% load data
allcfiles = list.files(path = "../data/countsf/", pattern = "Y(.+)counts$")


all_dfs = lapply(X = allcfiles, FUN = function(ff) {
  df = read.table(file = paste0("../data/countsf/", ff), header = T, sep = "\t", stringsAsFactors = F)
  
  cheader = strsplit(x = names(df)[7], split = "\\.")[[1]]
  cheader = paste0(strsplit(x = cheader[6], split = "_")[[1]][1:2], collapse = "_")

  names(df)[7] <- cheader
  
  df <- df[,c(1,7)]
  
  names(df)[1] <- "ensembl_gene_id"
  
  df$ensembl_gene_id <- sapply(X = strsplit(x = df$ensembl_gene_id, split = ":"), FUN = function(x) x[2])
  df <- as.data.table(df)
  setkey(df, "ensembl_gene_id")
  df
})

all_counts <- all_dfs[[1]]

for (ii in c(2:length(all_dfs))) {
  all_counts <- merge.data.table(x = all_counts, y = all_dfs[[ii]])
}


#%% make coldata
coldata = data.frame(sample=names(all_counts)[c(2:ncol(all_counts))], stringsAsFactors = F)

coldata$condition_num <- as.numeric(stringr::str_extract(coldata$sample, "\\d+(?=_)"))

# Extract second number (1, 2, or 3)
coldata$replicate_num <- stringr::str_extract(coldata$sample, "(?<=_)\\d+")

coldata$genotype <- dplyr::case_when(
  between(coldata$condition_num, 8, 11) ~ "WT",
  between(coldata$condition_num, 12, 15) ~ "Notch1KO",
  between(coldata$condition_num, 16, 19) ~ "p53KO",
  TRUE ~ NA_character_
)

# sort by condition_num, and then replicate_num
coldata <- coldata[order(coldata$condition_num, coldata$replicate_num), ]

coldata$treatment <- c('0' = 'Control', '1' = 'Hypoxia', '2' = 'Lactate', '3' = 'Oscillation')[ as.character(coldata$condition_num %% 4)]

coldata$condition <- paste0(coldata$genotype, '_', coldata$treatment )

# specify Control as the reference in treatment
coldata$treatment <- factor(coldata$treatment, levels = c("Control", "Hypoxia", "Lactate", "Oscillation"))

rownames(coldata) <- coldata$sample

#%% make count matrix file
count_mat <- as.matrix(x = all_counts[, c(2:ncol(all_counts)), with=F])
rownames(count_mat) <- all_counts$ensembl_gene_id

#%%

# re-order the columns of the count matrix to match the coldata
count_mat <- count_mat[, coldata$sample]

dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = coldata, design = ~ condition)

dds <- DESeq(dds)

#%% get PCA
vsd <- varianceStabilizingTransformation(dds)
pcaData <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$condition <- factor(pcaData$condition, levels = c("WT_Control", "WT_Hypoxia", "WT_Lactate", "WT_Oscillation",
                                                            "p53KO_Control", "p53KO_Hypoxia", "p53KO_Lactate", "p53KO_Oscillation",
                                                            "Notch1KO_Control", "Notch1KO_Hypoxia", "Notch1KO_Lactate", "Notch1KO_Oscillation"))
tiff("../figures/PCA_overall.tiff", height = 15, width = 15, units = "in", res = 300, compress = "lzw+p")
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = condition), size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33", "#A65628", "#999999", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#D95F02"))
dev.off()


#%%
# read in lengths from the ../data/countsf/YL9_1_counts.txt

lengths <- read.table(file = "../data/countsf/YL9_1.counts", header = T, sep = "\t", stringsAsFactors = F)
lengths$ensembl_gene_id <- sapply(X = strsplit(x = lengths$Geneid, split = ":"), FUN = function(x) x[2])

# create a data.table with gene IDs and lengths
lengths <- data.table(ensembl_gene_id = lengths$ensembl_gene_id, 
                     Length = lengths$Length)
setkey(lengths, "ensembl_gene_id")

# Calculate TPM
# First, create a function to calculate TPM
calculate_tpm <- function(counts, lengths) {
  # Step 1: Divide counts by length to get reads per kilobase (RPK)
  rpk <- counts / (lengths / 1000)
  scaling_factors <- colSums(rpk) / 1e6
  tpm <- sweep(rpk, 2, scaling_factors, "/")
  return(tpm)
}

# Calculate TPM for each sample
tpm_mat <- calculate_tpm(count_mat, lengths$Length)

# Convert to data.frame if needed
tpm_df <- as.data.frame(tpm_mat)
tpm_df$ensembl_gene_id <- rownames(tpm_df)

# merge with the namedf
tpm_df <- merge(tpm_df, namedf, by = "ensembl_gene_id", all.x = T)

# put the hgnc_symbol column first
tpm_df <- tpm_df[, c("hgnc_symbol", setdiff(colnames(tpm_df), "hgnc_symbol"))]

# convert all the othr column names except hgnc_symbol and ensembl_gene_id to corresponding condition_replicate from coldata
colnames(tpm_df)[c(-1,-2)] <- sapply(X = colnames(tpm_df)[c(-1,-2)], FUN = function(x) {
  cond = match(x, coldata$sample)
  paste0(coldata$condition[cond], "_", coldata$replicate_num[cond])
})

# write to an excel file
openxlsx::write.xlsx(tpm_df, file = "../results/tpm_mat.xlsx", rowNames = F, asTable = T, bandedRows = T)

#%% make separate dds for each genotype

dds_WT <- DESeqDataSetFromMatrix(countData = count_mat[, coldata$genotype == "WT"], colData = coldata[coldata$genotype == "WT", ], design = ~ treatment)
dds_WT <- DESeq(dds_WT)
vsd_WT <- varianceStabilizingTransformation(dds_WT)
pcaData_WT <- plotPCA(vsd_WT, intgroup = c("condition"), returnData = TRUE)
percentVar_WT <- round(100 * attr(pcaData_WT, "percentVar"))
pcaData_WT$condition <- factor(pcaData_WT$condition, levels = c("WT_Control", "WT_Hypoxia", "WT_Lactate", "WT_Oscillation"))
tiff("../figures/PCA_WT.tiff", height = 5, width = 5, units = "in", res = 300, compress = "lzw+p")
ggplot(pcaData_WT, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = name), size = 3) +
  xlab(paste0("PC1: ", percentVar_WT[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_WT[2], "% variance")) +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00"))
dev.off()

dds_p53KO <- DESeqDataSetFromMatrix(countData = count_mat[, coldata$genotype == "p53KO"], colData = coldata[coldata$genotype == "p53KO", ], design = ~ treatment)
dds_p53KO <- DESeq(dds_p53KO)
vsd_p53KO <- varianceStabilizingTransformation(dds_p53KO)
pcaData_p53KO <- plotPCA(vsd_p53KO, intgroup = c("condition"), returnData = TRUE)
percentVar_p53KO <- round(100 * attr(pcaData_p53KO, "percentVar"))
pcaData_p53KO$condition <- factor(pcaData_p53KO$condition, levels = c("p53KO_Control", "p53KO_Hypoxia", "p53KO_Lactate", "p53KO_Oscillation"))
tiff("../figures/PCA_p53KO.tiff", height = 5, width = 5, units = "in", res = 300, compress = "lzw+p")
ggplot(pcaData_p53KO, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = name), size = 3) +
  xlab(paste0("PC1: ", percentVar_p53KO[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_p53KO[2], "% variance")) +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00"))
dev.off()

dds_Notch1KO <- DESeqDataSetFromMatrix(countData = count_mat[, coldata$genotype == "Notch1KO"], colData = coldata[coldata$genotype == "Notch1KO", ], design = ~ treatment)
dds_Notch1KO <- DESeq(dds_Notch1KO)
vsd_Notch1KO <- varianceStabilizingTransformation(dds_Notch1KO)
pcaData_Notch1KO <- plotPCA(vsd_Notch1KO, intgroup = c("condition"), returnData = TRUE)
percentVar_Notch1KO <- round(100 * attr(pcaData_Notch1KO, "percentVar"))
pcaData_Notch1KO$condition <- factor(pcaData_Notch1KO$condition, levels = c("Notch1KO_Control", "Notch1KO_Hypoxia", "Notch1KO_Lactate", "Notch1KO_Oscillation"))
tiff("../figures/PCA_Notch1KO.tiff", height = 5, width = 5, units = "in", res = 300, compress = "lzw+p")
ggplot(pcaData_Notch1KO, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = name), size = 3) +
  xlab(paste0("PC1: ", percentVar_Notch1KO[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_Notch1KO[2], "% variance")) +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00"))
dev.off()


#%% get namedf df ensembl_gene_id to hgnc_symbol from biomart

if(! file.exists("namedf.rds")) {
  library(biomaRt)
  
  # Function to try different mirrors
  try_mirrors <- function() {
    mirrors <- c("useast", "uswest", "asia")
    for (m in mirrors) {
      tryCatch({
        mart <- useEnsembl(biomart = "ensembl", 
                          dataset = "hsapiens_gene_ensembl", 
                          mirror = m)
        message(paste("Successfully connected to", m, "mirror"))
        return(mart)
      }, error = function(e) {
        message(paste("Failed to connect to", m, "mirror"))
      })
    }
    stop("Failed to connect to any mirror")
  }
  
  mart <- try_mirrors()
  
  namedf <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  mart = mart)
  namedf <- as.data.table(namedf)
  setkey(namedf, "ensembl_gene_id")
  saveRDS(namedf, file = "../data/namedf.rds")
} else {
  namedf <- readRDS("../data/namedf.rds")
}


#%% compile the DE results and combine into one file each for each genotype, with columns like pvalue_condition, padj_condition, log2FoldChange_condition, stat_condition
# for each genotype

condition_comparisons <- list(c("treatment", "Hypoxia", "Control"),
                              c("treatment", "Lactate", "Control"),
                              c("treatment", "Oscillation", "Control"),
                              c("treatment", "Lactate", "Hypoxia"),
                              c("treatment", "Oscillation", "Hypoxia"),
                              c("treatment", "Oscillation", "Lactate"))

# for each genotype
for (genotype in c("WT", "p53KO", "Notch1KO")) {
  # get the dds for the genotype
  dds <- get(paste0("dds_", genotype))
  
  # get the results for each comparison
  res_list <- lapply(condition_comparisons, function(comp) {
    res <- results(dds, contrast = comp)
    suffix <- paste0(comp[2], "_vs_", comp[3])
    res <- as.data.frame(res)
    res$ensembl_gene_id <- rownames(res)

    # Keep only specific columns
    cols_to_keep <- c("ensembl_gene_id", 
                      if(identical(comp, condition_comparisons[[1]])) "baseMean" else NULL,  # Keep baseMean only from first comparison
                      "padj", "pvalue", "log2FoldChange", "stat")
    res <- res[, cols_to_keep]
    
    # add suffix to padj, pvalue, log2FoldChange, stat
    cols_to_rename <- intersect(colnames(res), c("padj", "pvalue", "log2FoldChange", "stat"))
    colnames(res)[colnames(res) %in% cols_to_rename] <- paste0(cols_to_rename, "_", suffix)
    res
  })
  
  # combine the results into one data frame
  res_combined <- Reduce(function(x, y) merge(x, y, by = "ensembl_gene_id"), res_list)
  
  # rename the columns
  colnames(res_combined)[-1] <- paste0(colnames(res_combined)[-1], "_", genotype)
  
  # merge with the namedf
  res_combined <- merge(res_combined, namedf, by = "ensembl_gene_id", all.x = T)
  
  # put the hgnc_symbol column first
  res_combined <- res_combined[, c("hgnc_symbol", setdiff(colnames(res_combined), "hgnc_symbol"))]


  # write the results to a xlsx file
  openxlsx::write.xlsx(res_combined, file = paste0("../results/DESeq2_results_", genotype, ".xlsx"), rowNames = F, asTable = T)
}

#%% do normoxia only DE analysis, p53KO vs WT, and Notch1KO vs WT

cdata_Normoxia <- coldata[coldata$treatment == "Control", ]
# specify WT as the reference in genotype
cdata_Normoxia$genotype <- factor(cdata_Normoxia$genotype, levels = c("WT", "p53KO", "Notch1KO"))

dds_Normoxia <- DESeqDataSetFromMatrix(countData = count_mat[, cdata_Normoxia$sample], colData = cdata_Normoxia, design = ~ genotype)
dds_Normoxia <- DESeq(dds_Normoxia)
vsd_Normoxia <- varianceStabilizingTransformation(dds_Normoxia)
pcaData_Normoxia <- plotPCA(vsd_Normoxia, intgroup = c("genotype"), returnData = TRUE)
percentVar_Normoxia <- round(100 * attr(pcaData_Normoxia, "percentVar"))
pcaData_Normoxia$genotype <- factor(pcaData_Normoxia$genotype, levels = c("WT", "p53KO", "Notch1KO"))

tiff("../figures/PCA_Normoxia.tiff", height = 5, width = 5, units = "in", res = 300, compress = "lzw+p")
ggplot(pcaData_Normoxia, aes(x = PC1, y = PC2, color = genotype)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = name), size = 3) +
  xlab(paste0("PC1: ", percentVar_Normoxia[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_Normoxia[2], "% variance")) +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"))
dev.off()


#%% write the DE results for normoxia to a xlsx file

numerator_list <- c("p53KO", "Notch1KO")
denominator <- "WT"

# get the results for each numerator, merge together, merge with namedf, and write to a xlsx file

res_list <- lapply(numerator_list, function(numerator) {
  res <- results(dds_Normoxia, contrast = c("genotype", numerator, denominator))
  res <- as.data.frame(res)
  res$ensembl_gene_id <- rownames(res)
  res
})

res_combined <- merge(res_list[[1]], res_list[[2]], by = "ensembl_gene_id", suffixes = paste0(".", numerator_list, "_vs_", denominator))

# merge with the namedf
res_combined <- merge(res_combined, namedf, by = "ensembl_gene_id", all.x = T)

# put the hgnc_symbol column first
res_combined <- res_combined[, c("hgnc_symbol", setdiff(colnames(res_combined), "hgnc_symbol"))]

head(res_combined)

# merge the namedf with the res_combined
res_combined <- merge(res_combined, namedf, by = "ensembl_gene_id", all.x = T)

# put the hgnc_symbol column first
res_combined <- res_combined[, c("hgnc_symbol", setdiff(colnames(res_combined), "hgnc_symbol"))]

# write to a xlsx file
openxlsx::write.xlsx(res_combined, file = "../results/DESeq2_results_Normoxia.xlsx", rowNames = F, asTable = T)








