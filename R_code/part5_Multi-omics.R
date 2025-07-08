################ Multi-omics data correlation with skin phenotypes ##########
# install.packages("devtools")
# devtools::install_github("hannet91/ggcor")

#1. Create files: Based on significantly changed species (17), significantly changed genes (350), significantly changed metabolites (3483)
file_names <- c("results/multiomics_analysis/mOTUs3_species_with_signif.txt", 
                "results/multiomics_analysis/kegg_genes_with_signif.txt", 
                "results/multiomics_analysis/metabolites_with_signif.txt")
for (file in file_names) {
  data <- read.csv(file, header = TRUE, check.names = FALSE, sep = "\t")
  transposed_data <- t(data)
  original_file_name <- basename(file)
  new_file_name <- paste0("transposed_", original_file_name)
  write.table(transposed_data, file.path("results/multiomics_analysis", new_file_name), sep = "\t", col.names = NA, quote = FALSE)
}
#Manually merge into: all_with_signif_data_for_mantel_gui_new.txt

#2. Read files
# table <- read.csv("all_data_for_mantel.txt", header=T, row.names =1,  sep="\t", check.names = F)#All microbes, genes, metabolites
table <- read.csv("results/multiomics_analysis/all_with_signif_data_for_mantel_gui_new.txt", header=T, row.names =1,  sep="\t", check.names = F)#Significantly changed microbes, genes, metabolites
# table <- read.csv("all_with_signif_data_for_mantel_gui2.txt", header=T, row.names =1,  sep="\t", check.names = F)#Not using day0 as baseline, microbes, genes, metabolites with at least 3 significant changes
# table <- read.csv("all_with_signif_data_for_mantel.txt", header=T, row.names =1,  sep="\t", check.names = F)
table<-t(table)
dim(table)

#skin_traits
skin_trait_table <- read.csv("skin_trait.csv", header=T, row.names =1, sep=",", check.names = F) 
View(skin_trait_table)
row.names(skin_trait_table)
dim(skin_trait_table)
head(skin_trait_table)
env_All <- skin_trait_table[-c(27),-c(1:2)]
dim(env_All)
row.names(env_All)
identical(rownames(table), rownames(env_All))
head(env_All)

#Match row names of the two matrices
env_All <- env_All[rownames(table) %in% rownames(env_All),]

#Environmental factor standardization
env_All_scale=scale(env_All)

#Calculate mantel and plot results
mantel_all <- fortify_mantel(table, env_All_scale, 
                             spec.select = list(Species = 1:17, Genes = 18:367, Metabolites = 368:3850), na.rm = T)
# mantel_all <- fortify_mantel(table, env_All_scale, 
#                              spec.select = list(Species = 1:13, Genes = 14:243, Metabolites = 244:1211), na.rm = T)
# mantel_all <- fortify_mantel(table, env_All_scale, 
#                              spec.select = list(Species = 1:13, Metabolites = 209:869), na.rm = T) 
print(mantel_all, n = 100)

mantel_all1 <- mantel_all %>% 
  mutate(r = cut(r, breaks = c(-Inf, 0.1, 0.2, Inf), 
                 labels = c("<0.1", "0.1-0.2", ">=0.2"),
                 right = FALSE),
         p.value = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                       labels = c("< 0.01", "0.01 - 0.051", ">= 0.051"),
                       right = FALSE))
print(mantel_all1, n = 100)

mantel_plot_all <- quickcor(env_All_scale, type = "upper", show.diag = T) + 
  geom_square() + 
  add_link(mantel_all1, mapping = aes(colour = p.value, size = r)) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c( "firebrick", "royalblue", "grey")) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))
mantel_plot_all  # Save as: mantel_plot_all

write.csv(mantel_all, "results/multiomics_analysis/all_with_signif_mantel_results.csv")

## Correlation heatmap of significantly changed genes, species with significantly changed metabolites
# 1. Genes-metabolites
# 2. Species-metabolites

# Create files: significantly changed species (17) + persistently significantly changed genes (5) + metabolites of interest (8), Correlation heatmap data
    # Create metabolite file enriched and baseline lower than tolerance group
    metab <- read.csv("metabolites.csv", header = TRUE, check.names = FALSE)
    metabolites_to_extract <- c(
      "Sphinganine",
      "SM(d18:0/20:2(11Z,14Z))",
      "Malic Acid",
      "Aldosterone 18-Glucuronide",
      "PC(22:4(7Z,10Z,13Z,16Z)/22:5(4Z,7Z,10Z,13Z,16Z))",
      "Dihydroxyacetone Phosphate Acyl Ester",
      "PC(22:6(4Z,7Z,10Z,13Z,16Z,19Z)/0:0)",
      "(9s,10s)-9,10-Dihydroxyoctadecanoic Acid"
    )
    selected_data <- metab[, c("sample", "Group", "Volunteer_ID", metabolites_to_extract)]
    order_vector <- c(
      # Baseline group
      "d0_01", "d0_03", "d0_06", "d0_07", "d0_11", "d0_15", "d0_23", "d0_27", "d0_35",
      # Adverse_reaction_occured group
      "d7_01", "d7_03", "d7_06", "d7_07", "d7_11", "d7_15", "d7_23", "d7_27", "d7_35",
      # Tolerance_established group (note order as given)
      "d21_01", "d21_03", "d28_06", "d28_07", "d28_11", "d21_15", "d21_23", "d28_27"
    )
    selected_data$sample <- factor(selected_data$sample, levels = order_vector)
    selected_data <- selected_data[order(selected_data$sample), ]
    selected_data$sample <- as.character(selected_data$sample)
    write.table(selected_data, file = "results/multiomics_analysis/specified_metabolites.txt", sep = "\t", row.names = FALSE, quote = FALSE)

  file_names <- c("results/multiomics_analysis/mOTUs3_species_with_signif.txt", 
                  "results/multiomics_analysis/kegg_genes_with_signify_both_2_stage.txt", 
                  "results/multiomics_analysis/specified_metabolites.txt")
  for (file in file_names) {
    data <- read.csv(file, header = TRUE, check.names = FALSE, sep = "\t")
    transposed_data <- t(data)
    original_file_name <- basename(file)
    new_file_name <- paste0("transposed_", original_file_name)
    write.table(transposed_data, file.path("results/multiomics_analysis", new_file_name), sep = "\t", col.names = NA, quote = FALSE)
  }
  #Manually merge into: Correlation_heatmap_data.txt


# Read data
df_raw <- read.csv("results/multiomics_analysis/Correlation_heatmap_data.txt",
                   header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
dim(df_raw)
head(df_raw[,1:5])  # View first few columns

# Divide data
df_species    <- df_raw[1:17, , drop = FALSE]     # Species
df_genes      <- df_raw[18:22, , drop = FALSE]    # Genes
df_metabolite <- df_raw[23:30, , drop = FALSE]    # Metabolites

# Function to calculate correlation coefficient matrix (returns r values)
calc_cor_matrix <- function(mat_x, mat_y, method = "spearman") {
  nX <- nrow(mat_x)
  nY <- nrow(mat_y)
  cor_mat <- matrix(NA, nrow = nX, ncol = nY)
  rownames(cor_mat) <- rownames(mat_x)
  colnames(cor_mat) <- rownames(mat_y)
  for(i in seq_len(nX)) {
    for(j in seq_len(nY)) {
      x <- as.numeric(mat_x[i, ])
      y <- as.numeric(mat_y[j, ])
      # Calculate correlation coefficient (pairwise.complete.obs removes NA)
      cor_val <- cor(x, y, method = method, use = "pairwise.complete.obs")
      cor_mat[i, j] <- cor_val
    }
  }
  return(cor_mat)
}

# Function to calculate p-value matrix
calc_pval_matrix <- function(mat_x, mat_y, method = "spearman") {
  nX <- nrow(mat_x)
  nY <- nrow(mat_y)
  pval_mat <- matrix(NA, nrow = nX, ncol = nY)
  rownames(pval_mat) <- rownames(mat_x)
  colnames(pval_mat) <- rownames(mat_y)
  for(i in seq_len(nX)) {
    for(j in seq_len(nY)) {
      x <- as.numeric(mat_x[i, ])
      y <- as.numeric(mat_y[j, ])
      # Use cor.test to calculate correlation coefficient and p-value
      test <- cor.test(x, y, method = method, use = "pairwise.complete.obs")
      pval_mat[i, j] <- test$p.value
    }
  }
  return(pval_mat)
}

# Convert to numeric matrices (ensure data is numeric)
mat_species <- as.matrix(df_species)
mat_metab   <- as.matrix(df_metabolite)
mat_genes   <- as.matrix(df_genes)

# Create significance marker function
create_sig_matrix <- function(pval_mat) {
  sig_mat <- matrix("", nrow = nrow(pval_mat), ncol = ncol(pval_mat))
  rownames(sig_mat) <- rownames(pval_mat)
  colnames(sig_mat) <- colnames(pval_mat)
  sig_mat[pval_mat < 0.05] <- "*"
  sig_mat[pval_mat < 0.01] <- "**"
  sig_mat[pval_mat < 0.001] <- "***"
  return(sig_mat)
}



#### 1. Species vs Metabolites ####

# Calculate correlation r value matrix
cor_mat_species_metab <- calc_cor_matrix(mat_species, mat_metab, method = "spearman")
# Calculate p-value matrix
pval_mat_species_metab <- calc_pval_matrix(mat_species, mat_metab, method = "spearman")

# Process matrix: remove all-NA rows/columns (if any), fill remaining NAs with 0
remove_rows <- which(rowSums(is.na(cor_mat_species_metab)) == ncol(cor_mat_species_metab))
if (length(remove_rows) > 0) {
  cor_mat_species_metab <- cor_mat_species_metab[-remove_rows, , drop=FALSE]
  pval_mat_species_metab <- pval_mat_species_metab[-remove_rows, , drop=FALSE]
}
remove_cols <- which(colSums(is.na(cor_mat_species_metab)) == nrow(cor_mat_species_metab))
if (length(remove_cols) > 0) {
  cor_mat_species_metab <- cor_mat_species_metab[, -remove_cols, drop=FALSE]
  pval_mat_species_metab <- pval_mat_species_metab[, -remove_cols, drop=FALSE]
}
cor_mat_species_metab[is.na(cor_mat_species_metab)] <- 0
pval_mat_species_metab[is.na(pval_mat_species_metab)]   <- 0

# Create significance marker matrix
sig_mat_species_metab <- create_sig_matrix(pval_mat_species_metab)

# Plot species-metabolite correlation heatmap (r values) with significance markers
pheatmap(
  cor_mat_species_metab,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Correlation: Species vs Metabolites (Spearman)",
  border_color = NA,
  display_numbers = sig_mat_species_metab,
  number_color = "black",
  fontsize_number = 10
) # Species-metabolite correlation heatmap

# Save correlation matrix as TSV
write.table(cor_mat_species_metab, file = "results/multiomics_analysis/species_metabolites_correlation.tsv", 
            sep = "\t", quote = FALSE, row.names = TRUE)
# Save p-value matrix as TSV
write.table(pval_mat_species_metab, file = "results/multiomics_analysis/species_metabolites_correlation_pvalues.tsv", 
            sep = "\t", quote = FALSE, row.names = TRUE)



#### 2. Genes vs Metabolites ####

# Calculate correlation r value matrix (genes vs metabolites)
cor_mat_gene_metab <- calc_cor_matrix(mat_genes, mat_metab, method = "spearman")
# Calculate p-value matrix
pval_mat_gene_metab <- calc_pval_matrix(mat_genes, mat_metab, method = "spearman")

# Remove all-NA rows/columns and fill NAs with 0
remove_rows <- which(rowSums(is.na(cor_mat_gene_metab)) == ncol(cor_mat_gene_metab))
if (length(remove_rows) > 0) {
  cor_mat_gene_metab <- cor_mat_gene_metab[-remove_rows, , drop=FALSE]
  pval_mat_gene_metab <- pval_mat_gene_metab[-remove_rows, , drop=FALSE]
}
remove_cols <- which(colSums(is.na(cor_mat_gene_metab)) == nrow(cor_mat_gene_metab))
if (length(remove_cols) > 0) {
  cor_mat_gene_metab <- cor_mat_gene_metab[, -remove_cols, drop=FALSE]
  pval_mat_gene_metab <- pval_mat_gene_metab[, -remove_cols, drop=FALSE]
}
cor_mat_gene_metab[is.na(cor_mat_gene_metab)] <- 0
pval_mat_gene_metab[is.na(pval_mat_gene_metab)]   <- 0

# Create significance marker matrix
sig_mat_gene_metab <- create_sig_matrix(pval_mat_gene_metab)

# Plot gene-metabolite correlation heatmap (r values) with significance markers
pheatmap(
  cor_mat_gene_metab,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Correlation: Genes vs Metabolites (Spearman)",
  border_color = NA,
  display_numbers = sig_mat_gene_metab,
  number_color = "black",
  fontsize_number = 10
) # Gene-metabolite correlation heatmap

# Save correlation matrix
write.table(cor_mat_gene_metab, file = "results/multiomics_analysis/genes_metabolites_correlation.tsv", 
          sep = "\t", quote = FALSE, row.names = TRUE)
# Save p-value matrix
write.table(pval_mat_gene_metab, file = "results/multiomics_analysis/genes_metabolites_correlation_pvalues.tsv", 
          sep = "\t", quote = FALSE, row.names = TRUE)