###############    KEGG_gene analysis     ###########################
#Gene diversity
species <- read.csv("kegg_genes.txt", header =T, check.names = F, sep = "\t")
dim(species)
head(species)
species$sample <- factor(species$sample)
species$Group <- factor(species$Group)
species$Volunteer_ID <- factor(species$Volunteer_ID)
colnames(species[,-c(1:3)])[1]

#shannon
shannon <-as.matrix(diversity(species[,-c(1:3)],index="shannon"))
shannon
shannon_table <- cbind.data.frame(species[,1:3], shannon)
shannon_table
write.csv(shannon_table, "results/gene_analysis/Genes_shannon_index_table.csv")
dim(shannon_table)
head(shannon_table)
sink("results/gene_analysis/Genes_shannon_results.csv", split = T)
#calculate mean_sd
shannon_results <- shannon_table %>% group_by(Group) %>% summarise(mean = mean(shannon), sd = sd(shannon))
print(shannon_results)
#calculate significance
# Remove rows where sample column is d0_35 and d7_35
shannon_table <- shannon_table %>% filter(!(sample %in% c("d0_35", "d7_35")))
shannon_stat_test <- compare_means(shannon ~ Group, 
                                   paired = T, 
                                   data = shannon_table,
                                   p.adjust.method = "BH")
print(shannon_stat_test)
sink()

#richness
richness <-as.matrix(specnumber(species[,-c(1:3)]))
richness
richness_table <- cbind.data.frame(species[,1:3], richness)
richness_table
write.csv(richness_table, "results/gene_analysis/Genes_richness_table.csv") 
dim(richness_table)
head(richness_table)
sink("results/gene_analysis/Genes_richness_results.csv", split = T)
#calculate mean_sd
richness_results <- richness_table %>% group_by(Group) %>% summarise(mean = mean(richness), sd = sd(richness))
print(richness_results)
#calculate significance
# Remove rows where sample column is d0_35 and d7_35
richness_table <- richness_table %>% filter(!(sample %in% c("d0_35", "d7_35")))
richness_stat_test <- compare_means(richness ~ Group, 
                                    paired = T, 
                                    data = richness_table,
                                    p.adjust.method = "BH")
print(richness_stat_test)
sink()

#evenness
evenness <- as.matrix(diversity(species[,-c(1:3)],index="shannon")/log(specnumber(species[,-c(1:3)])))
evenness
evenness_table <- cbind.data.frame(species[,1:3], evenness)
evenness_table
write.csv(evenness_table, "results/gene_analysis/Genes_evenness_table.csv") 
dim(evenness_table)
head(evenness_table)
sink("results/gene_analysis/Genes_evenness_results.csv", split = T)
#calculate mean_sd
evenness_results <- evenness_table %>% group_by(Group) %>% summarise(mean = mean(evenness), sd = sd(evenness))
print(evenness_results)
#calculate significance
# Remove rows where sample column is d0_35 and d7_35
evenness_table <- evenness_table %>% filter(!(sample %in% c("d0_35", "d7_35")))
evenness_stat_test <- compare_means(evenness ~ Group, 
                                    paired = T, 
                                    data = evenness_table,
                                    p.adjust.method = "BH")
print(evenness_stat_test)
sink()

#plotting
comparisons <- list(c("Baseline", "Adverse_reaction_occured"), 
                    c("Baseline", "Tolerance_established"),
                    c("Adverse_reaction_occured", "Tolerance_established"))
# Shannon violin plot
shannon_violinplot <- ggviolin(shannon_table, x = "Group", y = "shannon",
                               color = "Group", fill = "Group", palette = c("grey", "#c97937", "#3b6291"), add = "mean_sd",
                               add.params = list(color = "white", size = 0.5)) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif", paired = TRUE, p.adjust.method = "BH") +
  labs(y = "Species Shannon index", x = "") +
  theme_bw()
shannon_violinplot
# Shannon line plot
shannon_lineplot <- ggline(shannon_table, x = "Group", y = "shannon", group = "Volunteer_ID",
                           size = 0.3, color = "Volunteer_ID", palette = "lancet") +
  labs(x = NULL, y = "Species Shannon index") +
  theme_bw()
shannon_lineplot
# Richness violin plot
richness_violinplot <- ggviolin(richness_table, x = "Group", y = "richness",
                                color = "Group", fill = "Group", palette = c("grey", "#c97937", "#3b6291"), add = "mean_sd",
                                add.params = list(color = "white", size = 0.5)) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif", paired = TRUE, p.adjust.method = "BH") +
  labs(y = "Species richness", x = "") +
  theme_bw()
# Richness line plot
richness_lineplot <- ggline(richness_table, x = "Group", y = "richness", group = "Volunteer_ID",
                            size = 0.3, color = "Volunteer_ID", palette = "lancet") +
  labs(x = NULL, y = "Species richness") +
  theme_bw()
richness_lineplot
# Evenness violin plot
evenness_violinplot <- ggviolin(evenness_table, x = "Group", y = "evenness",
                                color = "Group", fill = "Group", palette = c("grey", "#c97937", "#3b6291"), add = "mean_sd",
                                add.params = list(color = "white", size = 0.5)) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif", paired = TRUE, p.adjust.method = "BH") +
  labs(y = "Species evenness", x = "") +
  theme_bw()
# Evenness line plot
evenness_lineplot <- ggline(evenness_table, x = "Group", y = "evenness", group = "Volunteer_ID",
                            size = 0.3, color = "Volunteer_ID", palette = "lancet") +
  labs(x = NULL, y = "Species evenness") +
  theme_bw()
evenness_lineplot
# Combined plots
alpha_diversity_violinplot <- ggarrange(
  ggarrange(shannon_violinplot, richness_violinplot, evenness_violinplot, 
            nrow = 3, ncol = 1, common.legend = T, legend = "right"),
  ggarrange(shannon_lineplot, richness_lineplot, evenness_lineplot, 
            nrow = 3, ncol = 1, common.legend = T, legend = "right"),
  nrow = 1, ncol = 2)
alpha_diversity_violinplot  # Save as: Gene α diversity

#PCoA plot for all genes
data = species
data[,c(1:3)]
data1 = species[,-c(1:3)]
colnames(data1)[1]
distmatrix <- vegdist(data1, method='bray')
pcoa<- pcoa(distmatrix, correction = "none", rn = NULL)
pcoa$vectors
dim(pcoa$vectors)
row.names(pcoa$vectors)
pcoa_new <- cbind.data.frame(data[,c(1:3)], pcoa$vectors)
head(pcoa_new)
#### Calculate PC1 and PC2
PC1 = pcoa$vectors[,1]
PC2 = pcoa$vectors[,2]
## Convert PC proportions
pc1 <-round(pcoa$values$Relative_eig[1]*100,2)
pc1
pc2 <-round(pcoa$values$Relative_eig[2]*100,2)
pc2
#Construct PCoA table containing PC1 and PC2 axes
species_pcoa <- pcoa_new[, c(2:5)]
head(species_pcoa)
colnames(species_pcoa)[3] <- "PC1"
colnames(species_pcoa)[4] <- "PC2"
dim(species_pcoa)
head(species_pcoa)
#plotting
species_pcoa$Group <- factor(species_pcoa$Group, levels = c("Baseline", "Adverse_reaction_occured", "Tolerance_established"), ordered = T)
species_pcoa_plot = ggscatter(species_pcoa, x = "PC1", y = "PC2", 
                                 color = "Group", palette = c("grey","#c97937","#3b6291"), ellipse = TRUE, ellipse.type = "norm",
                                 ellipse.border.remove = F, ellipse.alpha = 0.05, ellipse.level = 0.7, 
                                 size = 1,  legend = "right", mean.point =T, mean.point.size = 4,
                                 xlab = paste0("PCoA1 ( ",pc1,"%"," )",sep="" ),
                                 ylab = paste0("PCoA2 ( ",pc2,"%"," )",sep="" ))+theme_bw()
species_pcoa_plot   # Save as: Gene PCoA plot
######## adonis dissimilarity analysis
sink('results/gene_analysis/Genes_dissimilarity_test_results.txt', split = T)
set.seed(1234)
pairwise.adonis2(data1~Group, data, nperm = 9999)
sink()

#Significance testing for genes
mOTUs3_all_species <- read.csv("kegg_genes.txt", header = TRUE, check.names = FALSE, sep = "\t")
dim(mOTUs3_all_species)
head(mOTUs3_all_species)[1:6]
mOTUs3_all_species$sample <- factor(mOTUs3_all_species$sample)
mOTUs3_all_species$Volunteer_ID <- factor(mOTUs3_all_species$Volunteer_ID)
mOTUs3_all_species$Group <- factor(mOTUs3_all_species$Group)
colnames(mOTUs3_all_species[,-c(1:3)])[1]
species_melt_table <- melt(mOTUs3_all_species, id.vars = 1:3, variable.name = "Species", value.name = "Value")
species_melt_table$Value[species_melt_table$Value == 0] <- 0.00001  # Replace 0 values with a very small value
dim(species_melt_table)
head(species_melt_table, 20)
tail(species_melt_table, 20)
species_abundance_result <- species_melt_table %>% group_by(Species, Group) %>% summarise(mean = mean(Value), sd = sd(Value))
print(species_abundance_result, n = 100)
write.csv(species_abundance_result, "results/gene_analysis/Genes_abundance_mean_sd_Group.csv")
data1 <- species_melt_table
data1 <- data1 %>% filter(!(sample %in% c("d0_35", "d7_35")))
species_stat_test1 <- compare_means(Value ~ Group, group.by = c("Species"), paired = TRUE, data = data1, p.adjust.method = "BH")
print(species_stat_test1, n = 200)
write.csv(species_stat_test1, "results/gene_analysis/Genes_abundance_stat_test_Group.csv")

#Create significantly changed gene count and information file
# Step 1: Read file and extract significantly changed genes
gene_stat <- read.csv("results/gene_analysis/Genes_abundance_stat_test_Group.csv", header = TRUE, sep = ",")
# Step 2: Extract significantly changed genes where p.signif column is not "ns" and remove duplicates
significant_genes <- unique(gene_stat$Species[gene_stat$p.signif != "ns"])
# Step 3: Read kegg_genes.txt file
kegg_genes <- read.csv("kegg_genes.txt", header = TRUE, check.names = FALSE, sep = "\t")
# Step 4: Filter kegg_genes to keep only significantly changed gene columns
filtered_kegg_genes <- kegg_genes[, c("sample", "Group", "Volunteer_ID", significant_genes)]
# Step 5: Save filtered results to new file
write.table(filtered_kegg_genes, "kegg_genes_with_signif.txt", sep = "\t", row.names = FALSE, quote = FALSE)
head(filtered_kegg_genes)
# Step 6: Count significantly changed genes by comparison group and save gene information separately
gene_stat$comparison <- paste(gene_stat$group1, gene_stat$group2, sep = "_vs_")
# Count significantly changed genes by comparison group
comparisons <- unique(gene_stat$comparison)
for (comparison in comparisons) {
  comparison_genes <- unique(gene_stat$Species[gene_stat$comparison == comparison & gene_stat$p.signif != "ns"])
  write.csv(comparison_genes, paste0("results/gene_analysis/signif_genes_", comparison, ".csv"), row.names = FALSE)
}
# Step 7: Count significantly changed genes by group and save statistics
signif_summary <- gene_stat %>% 
  filter(p.signif != "ns") %>%
  group_by(group1, group2) %>%
  summarise(significant_gene_count = n(), genes = list(Species))
signif_summary$genes <- sapply(signif_summary$genes, function(x) paste(x, collapse = ";"))
write.csv(signif_summary, "results/gene_analysis/signif_gene_summary_by_group.csv", row.names = FALSE)


#There are 350 significantly changed genes, here we only plot genes that significantly changed in both phases
#1. Create significantly changed gene file
# Step 1: Read Genes_abundance_stat_test_Group.csv and extract significantly changed genes
gene_stat <- read.csv("results/gene_analysis/Genes_abundance_stat_test_Group.csv", header = TRUE, sep = ",")
signif_genes_baseline <- gene_stat$Species[gene_stat$group1 == "Adverse_reaction_occured" & 
                                           gene_stat$group2 == "Baseline" & 
                                           gene_stat$p.signif != "ns"]

signif_genes_tolerance <- gene_stat$Species[gene_stat$group1 == "Adverse_reaction_occured" & 
                                            gene_stat$group2 == "Tolerance_established" & 
                                            gene_stat$p.signif != "ns"]
# Extract genes that are significant in both comparison groups
signif_genes <- intersect(signif_genes_baseline, signif_genes_tolerance)
# Remove duplicate genes
signif_genes <- unique(signif_genes)
# Step 2: Read kegg_genes.txt file
kegg_genes <- read.csv("kegg_genes.txt", header = TRUE, sep = "\t", check.names = FALSE)
# Step 3: Filter kegg_genes to keep only significantly changed gene columns
filtered_kegg_genes <- kegg_genes[, c("sample", "Group", "Volunteer_ID", signif_genes)]
# Step 4: Save filtered results
write.table(filtered_kegg_genes, "kegg_genes_with_signify_both_2_stage.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#2. Plotting
skin_trait <- read.csv("kegg_genes_with_signify_both_2_stage.txt", header = TRUE, check.names = FALSE, sep = "\t")
skin_trait <- skin_trait[!(skin_trait$sample %in% c("d0_35", "d7_35")), ]
skin_trait$sample <- factor(skin_trait$sample)
skin_trait$Group <- factor(skin_trait$Group, levels = c("Baseline", "Adverse_reaction_occured", "Tolerance_established"), ordered = T)
skin_trait$Volunteer_ID <- factor(skin_trait$Volunteer_ID)
colnames(skin_trait[,-c(1:3)])[1]
#Group-wise statistical calculations
skin_trait_melt_table <- reshape2::melt(skin_trait, id = 1:3, variable.name = "skin_trait_index", value.name = "Value")
#plotting
comparisons <- list(c("Baseline", "Adverse_reaction_occured"), 
                    c("Adverse_reaction_occured", "Tolerance_established"),
                    c("Baseline", "Tolerance_established"))
skin_trait_violinplot <- ggboxplot(skin_trait_melt_table, x = "Group", y = "Value",
                                  color = "Group",fill = "Group",  palette = c("grey","#c97937","#3b6291"), add = "mean_sd",
                                  add.params = list(color = "black", size = 0.2))+
  theme_bw()+
  labs(x = "")+labs(y = "")+
  facet_wrap(~skin_trait_index, ncol = 1, scales = "free_y") +  # Set to display in one column
  theme(strip.text = element_text(size = 8)) +  # Adjust subtitle font size
  stat_compare_means(method = "wilcox.test", paired = TRUE, label = "p.signif", 
                    hide.ns = FALSE, label.y.npc = "top", vjust = 0.5,  comparisons = comparisons, p.adjust.method = "BH")

# Plot line graph with Volunteer_ID as x-axis
lineplot <- ggline(skin_trait_melt_table, x = "Group", y = "Value", 
                   group = "Volunteer_ID", color = "Volunteer_ID", palette = "lancet",
                   add = "mean_se", facet.by = "skin_trait_index", ncol = 1, scales = "free_y") +  # Right column, ncol = 1
  theme_bw() +
  labs(x = "Group", y = "") +  # x-axis is Group
  theme(strip.text = element_text(size = 8)) 

combined_plot <- ggarrange(skin_trait_violinplot, lineplot, ncol = 2)
combined_plot # Save as: Genes significantly changed in both phases


#PCoA plot for significantly changed genes
species <- read.csv("kegg_genes_with_signif.txt", header = TRUE, check.names = FALSE, sep = "\t")
row_sums <- rowSums(species)
species_filtered <- species[row_sums != 0, ]  # Filter all-zero rows
species <- na.omit(species_filtered)  # Filter missing values
data = species
data1 = species[,-c(1:3)]
colnames(data1)[1]
distmatrix <- vegdist(data1, method='bray')
pcoa<- pcoa(distmatrix, correction = "none", rn = NULL)
pcoa$vectors
dim(pcoa$vectors)
row.names(pcoa$vectors)
pcoa_new <- cbind.data.frame(data[,c(1:3)], pcoa$vectors)
head(pcoa_new)
#### Calculate PC1 and PC2
PC1 = pcoa$vectors[,1]
PC2 = pcoa$vectors[,2]
## Convert PC proportions
pc1 <-round(pcoa$values$Relative_eig[1]*100,2)
pc1
pc2 <-round(pcoa$values$Relative_eig[2]*100,2)
pc2
#Construct PCoA table containing PC1 and PC2 axes
species_pcoa <- pcoa_new[, c(2:5)]
head(species_pcoa)
colnames(species_pcoa)[3] <- "PC1"
colnames(species_pcoa)[4] <- "PC2"
dim(species_pcoa)
head(species_pcoa)
#plotting
species_pcoa$Group <- factor(species_pcoa$Group, levels = c("Baseline", "Adverse_reaction_occured", "Tolerance_established"), ordered = T)
species_pcoa_plot = ggscatter(species_pcoa, x = "PC1", y = "PC2", 
                                 color = "Group", palette = c("grey","#c97937","#3b6291"), ellipse = TRUE, ellipse.type = "norm",
                                 ellipse.border.remove = F, ellipse.alpha = 0.05, ellipse.level = 0.7, 
                                 size = 1,  legend = "right", mean.point =T, mean.point.size = 4,
                                 xlab = paste0("PCoA1 ( ",pc1,"%"," )",sep="" ),
                                 ylab = paste0("PCoA2 ( ",pc2,"%"," )",sep="" ))+theme_bw()
species_pcoa_plot   # Save as: PCoA plot of significantly changed genes
######## adonis dissimilarity analysis
sink('results/gene_analysis/significant_Genes_dissimilarity_test_results.txt', split = T)
set.seed(1234)
pairwise.adonis2(data1~Group, data, nperm = 9999)
sink()


####### Gene comparison analysis between 9 retinol-tolerant and 9 intolerant volunteers #######
# Inter-group gene significance differences
skin_trait <- read.csv("results/gene_analysis/kegg_gene_with_signif_retinol_yes_no_new.csv", header=T, sep=",", check.names = F) 
dim(skin_trait)
skin_trait$sample <- factor(skin_trait$sample)
skin_trait$Group <- factor(skin_trait$Group, levels = c("Intolerant_Baseline", "Tolerant_Baseline"), ordered = T)
colnames(skin_trait[,-c(1:2)])[1]
#Group-wise statistical calculations
skin_trait_melt_table <- reshape2::melt(skin_trait, id = 1:2, variable.name = "skin_trait_index", value.name = "Value")
dim(skin_trait_melt_table)
head(skin_trait_melt_table,18)
sink("results/gene_analysis/gene_with_signif_retinol_yes_no_results.csv", split = T)
#calculate mean_sd
skin_trait_mean_sd_results <- skin_trait_melt_table %>% group_by(skin_trait_index, Group) %>% summarise(mean = mean(Value), sd = sd(Value))
print(skin_trait_mean_sd_results, n=15000)
#calculate significance
skin_trait_stat_test <- compare_means(Value ~ Group, group.by = "skin_trait_index", 
                                      #paired = T,
                                      data = skin_trait_melt_table)
print(skin_trait_stat_test, n=10000)
sink()
#plotting
# skin_trait_melt_table$Value[skin_trait_melt_table$Value == 0] <- 0.00001 
selected_traits <- c("K00430: peroxidase [EC:1.11.1.7]", "K03949: NADH dehydrogenase (ubiquinone) 1 alpha subcomplex subunit 5", "K10365: capping protein (actin filament) muscle Z-line, beta", "K10850: MFS transporter, NNP family, putative nitrate transporter", "K11517: (S)-2-hydroxy-acid oxidase [EC:1.1.3.15]")  # Replace with actual feature names
#Plot boxplots for specified columns only
filtered_skin_trait <- skin_trait_melt_table %>%
  filter(skin_trait_index %in% selected_traits)
comparisons <- list(c("Intolerant_Baseline", "Tolerant_Baseline"))
skin_trait_violinplot <- ggboxplot(filtered_skin_trait, x = "Group", y = "Value",
                                  color = "Group",fill = "Group",  palette = c("#c97937","#3b6291"), add = "mean_sd",
                                  add.params = list(color = "black", size = 0.2))+
  theme_bw()+
  labs(x = "")+labs(y = "")+
  facet_wrap(~skin_trait_index, ncol = 1, nrow =5, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.8)) +
  stat_compare_means(method = "wilcox.test", paired = FALSE, label = "p.signif", 
                    hide.ns = FALSE, label.y.npc = "top", vjust = 0.5, comparisons = comparisons, p.adjust.method = "BH")
skin_trait_violinplot  # Save as: Inter-group baseline gene comparison


## Differential gene KEGG pathway analysis
#Prepare KO abundance table and experimental metadata

#Create retinol_KEGG_gene_for_reporter_score.csv
See get_kegg_report_data.py file

KO_abundance = read.csv("results/gene_analysis/retinol_KEGG_gene_for_reporter_score.csv", header = T, sep=',', check.names = F, row.names = 1)  #rows are KOs, columns are samples
metadata = read.csv("results/gene_analysis/retinol_group_for_reporter_score.csv", header = T, sep=',', check.names = F, row.names = 1)  #grouping information, first column is sample name, second column is group
load_KOlist()
head(KOlist$pathway)
metadata$Group <- factor(metadata$Group, levels = c("Tolerant_Baseline", "Intolerant_Baseline"), ordered = TRUE)
cat("Comparison: ",levels(factor(metadata$Group)))
#Note: write the compared group (control) at the end, whether the front group (treatment group) is upregulated or downregulated compared to the back group (control group)
#Because the R package has been updated: group as a factor variable, the first level will be set as the control group, you can change the factor level to change your comparison.
#In the code below, control is written at the end, so the reporter score calculation is reversed.
set.seed(1234)
# Change of the former (10:18 intolerant) relative to the latter (1:9 tolerant)
reporter_score_res1=reporter_score(KO_abundance[, c(10:18, 1:9)], metadata[c(10:18 ,1:9), ], mode="directed", method = "wilcox.test",
                                  p.adjust.method1 = "BH", type = "pathway")                             
View(reporter_score_res1$reporter_s)
dim(reporter_score_res1$reporter_s)
plot1 <- plot_report(reporter_score_res1,rs_threshold = 1.96, y_text_size = 5)
plot1   # retinol_Gene_enrichment_Baseline
reporter_score_res <- cbind(reporter_score_res1$reporter_s)
write.csv(reporter_score_res, "results/gene_analysis/retinol_reporter_score_results.csv") 

# Read data file
df <- read.csv("results/gene_analysis/retinol_reporter_score_results.csv")
filtered_df <- df[df$ReporterScore <= -1.96, ]
write.csv(filtered_df, "results/gene_analysis/retinol_RS_with_enrich.csv", row.names = FALSE)


#pcoa plot
#Environmental factor standardization
skin_trait <- read.csv("results/gene_analysis/kegg_gene_with_signif_retinol_yes_no.csv", header=T, sep=",", check.names = F) 
row_sums <- rowSums(skin_trait[, -c(1, 2)] )
skin_trait_filtered <- skin_trait[row_sums != 0, ]  # Filter all-zero rows
skin_trait <- na.omit(skin_trait_filtered)  # Filter missing values

data=skin_trait
data1=skin_trait[,-c(1:2)]
colnames(data1)[1]
distmatrix <- vegdist(data1, method='bray')
pcoa<- pcoa(distmatrix, correction = "none", rn = NULL)
pcoa$vectors
dim(pcoa$vectors)
row.names(pcoa$vectors)
pcoa_new <- cbind.data.frame(data[,c(1:2)], pcoa$vectors)
head(pcoa_new)
#### Calculate PC1 and PC2
PC1 = pcoa$vectors[,1]
PC2 = pcoa$vectors[,2]
## Convert PC proportions
pc1 <-round(pcoa$values$Relative_eig[1]*100,2)
pc1
pc2 <-round(pcoa$values$Relative_eig[2]*100,2)
pc2
#Construct PCoA table containing PC1 and PC2 axes
species_pcoa <- pcoa_new[, c(2:4)]
head(species_pcoa)
colnames(species_pcoa)[2] <- "PC1"
colnames(species_pcoa)[3] <- "PC2"
dim(species_pcoa)
head(species_pcoa)
#plotting
species_pcoa$Group <- factor(species_pcoa$Group, levels = c("Intolerant_Baseline", "Tolerant_Baseline"), ordered = T)
skin_trait_pcoa_plot = ggscatter(species_pcoa, x = "PC1", y = "PC2", 
                                 color = "Group", palette = c("#c97937","#3b6291"), ellipse = TRUE, ellipse.type = "norm",
                                 ellipse.border.remove = F, ellipse.alpha = 0.05, ellipse.level = 0.7, 
                                 size = 1,  legend = "right", mean.point =T, mean.point.size = 4,
                                 xlab = paste0("PCoA1 ( ",pc1,"%"," )",sep="" ),
                                 ylab = paste0("PCoA2 ( ",pc2,"%"," )",sep="" ))+theme_bw()
skin_trait_pcoa_plot  # Save as: Inter-group baseline gene comparison PCoA
#Dissimilarity analysis
sink('results/gene_analysis/gene_with_signif_retinol_yes_no_dissimilarity_test_results.txt', split = T)
set.seed(1234)
pairwise.adonis2(data1~Group, data, nperm = 9999)
sink()


############ KEGG enrichment analysis calculating reporter score ####################
#if(!require("devtools"))install.packages("devtools")
#devtools::install_github('Asa12138/pcutils',dependencies=T)
#devtools::install_github('Asa12138/ReporterScore',dependencies=T)
#Prepare KO abundance table and experimental metadata
KO_abundance = read.csv("results/gene_analysis/KEGG_gene_for_reporter_score.csv", header = T, sep=',', check.names = F, row.names = 1)  #rows are KOs, columns are samples
dim(KO_abundance)
head(KO_abundance[,1:6])
metadata = read.csv("results/gene_analysis/group_for_reporter_score.csv", header = T, sep=',', check.names = F, row.names = 1)  #grouping information, first column is sample name, second column is group
dim(metadata)
head(metadata)
load_KOlist()
head(KOlist$pathway)
cat("Comparison: ",levels(factor(metadata$Group)))
## Comparison:  OE WT
###KEGG gene analysis###
KO_pvalue=ko.test(KO_abundance,"Group", metadata, method = "kruskal.test", p.adjust.method = "fdr")
head(KO_pvalue)
KO_pvalue_L2FC <- KO_pvalue %>% mutate(L2FC_In_Ba = round(log2(average_Intolerance/average_Baseline),2), 
                                       L2FC_To_Ba = round(log2(average_Tolerance/average_Baseline),2),
                                       L2FC_To_In = round(log2(average_Tolerance/average_Intolerance),2))
head(KO_pvalue_L2FC)
write.csv(KO_pvalue_L2FC, "results/gene_analysis/KO_pvalue_L2FC_results.csv")
reporter_score_res1=reporter_score(KO_abundance[,1:18], "Group", metadata, mode="directed", method = "wilcox.test",
                                   p.adjust.method = "fdr", type = "pathway")
reporter_score_res2=reporter_score(KO_abundance[,-c(10:18)], "Group", metadata, mode="directed", method = "wilcox.test",
                                   p.adjust.method = "fdr", type = "pathway")
reporter_score_res3=reporter_score(KO_abundance[,10:26], "Group", metadata, mode="directed", method = "wilcox.test",
                                   p.adjust.method = "fdr", type = "pathway")
## ==================================load KOlist===================================
## ===================KOlist download time: 2023-05-12 00:07:41====================
## If you want to update KOlist, use `update_KO_file()`
## ===================================1.KO test====================================
## =========================2.Transfer p.value to z-score==========================
## ==========================3.Calculating reporter score==========================
View(reporter_score_res1$reporter_s)
reporter_score_res <- rbind.data.frame(reporter_score_res1$reporter_s, reporter_score_res2$reporter_s, reporter_score_res3$reporter_s)
write.csv(reporter_score_res1$reporter_s, "results/gene_analysis/reporter_score_results.csv")
plot <- plot_report(reporter_score_res,rs_threshold = 1.96, y_text_size = 5)
plot
#After manual modification and order adjustment
# Plot bar plot
KEGG_gene_RS <- read.csv("results/gene_analysis/KEGG_gene_RS.txt", header = T, sep='\t', check.names = F)
head(KEGG_gene_RS)
KEGG_gene_RS$Description <- factor(KEGG_gene_RS$Description,
                                   levels = as.vector(KEGG_gene_RS$Description),
                                   ordered = TRUE)  #Fix species order through levels and ordered
levels(factor(KEGG_gene_RS$Group))
p <- ggplot(KEGG_gene_RS, aes(x = ReporterScore, y = Description, fill = Group)) +
  geom_bar(stat = "identity") + xlim(-8, 10)+
  scale_fill_manual(values = c("#1F78B4", "#D2691E")) +
  labs(y = " ", x = paste0("Reporter Score")) +
  theme_classic()+
  geom_vline(xintercept = c(-1.96, 1.96), linetype = "dashed")+
  theme(legend.text = element_text(size = 8))
p  #KEGG_pathway_RS_plot

#When we focus on one pathway, e.g. "map00780":
plot_KOs_in_pathway(reporter_score_res,map_id = "map00780")
#display as a network:
## Loading required namespace: MetaNet
## Loading required namespace: tidyr
plot_KOs_network(reporter_score_res,map_id = "map00780",main="")
#look at the KOs abundance in a pathway:
## Loading required namespace: ggpubr
## `geom_smooth()` using formula = 'y ~ x'
plot_KOs_box(reporter_score_res,only_sig = TRUE)
#display as a heatmap:
## Loading required namespace: pheatmap
plot_KOs_heatmap(reporter_score_res,only_sig = TRUE,heatmap_param = list(cutree_rows=2))

#summaries each levels abundance
load_KO_htable()
View(KO_htable)
head(KO_htable)
plot_KO_htable()

KO_level_pathway=up_level_KO(KO_abundance,level = "pathway",show_name = TRUE)
pcutils::stackplot(KO_level_pathway[-which(rownames(KO_level_pathway)=="Unknown"),])

KO_level1=up_level_KO(KO_abundance,level = "level1",show_name = TRUE)
View(KO_level1)
pcutils::stackplot(KO_level1[-which(rownames(KO_level1)=="Unknown"),])
ko_pvalue=ko.test(KO_level1,"Group",metadata,p.adjust.method = "fdr")
ko_pvalue

KO_level2=up_level_KO(KO_abundance,level = "level2",show_name = TRUE)
pcutils::stackplot(KO_level2[-which(rownames(KO_level2)=="Unknown"),])


############ KEGG enrichment analysis calculating reporter score ####################
# if(!require("devtools"))install.packages("devtools")
# BiocManager::install("KEGGREST")
# devtools::install_github('Asa12138/pcutils',dependencies=T)
# devtools::install_github('Asa12138/ReporterScore',dependencies=T)

#Prepare KO abundance table and experimental metadata
KO_abundance = read.csv("results/gene_analysis/KEGG_gene_for_reporter_score.csv", header = T, sep=',', check.names = F, row.names = 1)  #rows are KOs, columns are samples
metadata = read.csv("results/gene_analysis/group_for_reporter_score.csv", header = T, sep=',', check.names = F, row.names = 1)  #grouping information, first column is sample name, second column is group
load_KOlist()
head(KOlist$pathway)

metadata$Group <- factor(metadata$Group, levels = c("Baseline", "Adverse_reaction_occured", "Tolerance_established"), ordered = TRUE)
cat("Comparison: ",levels(factor(metadata$Group)))

#Note: write the compared group (control) at the end, whether the front group (treatment group) is upregulated or downregulated compared to the back group (control group)
#Because the R package has been updated: group as a factor variable, the first level will be set as the control group, you can change the factor level to change your comparison.
#In the code below, control is written at the end, so the reporter score calculation is reversed.
set.seed(1234)
# Baseline vs Adverse_reaction_occured
# Change of the former (10:18 adverse reaction occurred) relative to the latter (1:9 baseline)
reporter_score_res1=reporter_score(KO_abundance[, c(10:18, 1:9)], metadata[c(10:18 ,1:9), ], mode="directed", method = "wilcox.test",
                                  p.adjust.method1 = "BH", type = "pathway")                             
View(reporter_score_res1$reporter_s)
dim(reporter_score_res1$reporter_s)
plot1 <- plot_report(reporter_score_res1,rs_threshold = 1.96, y_text_size = 5)
plot1   # Gene_enrichment_Baseline_Intolerance

set.seed(1234)
# Baseline vs Tolerance_established
# Change of the former (19:26 tolerance established) relative to the latter (1:9 baseline)
reporter_score_res2=reporter_score(KO_abundance[, c(19:26, 1:9)], metadata[c(19:26, 1:9), ], mode="directed", method = "wilcox.test",
                                   p.adjust.method1 = "BH", type = "pathway")
View(reporter_score_res2$reporter_s)
dim(reporter_score_res2$reporter_s)
plot2 <- plot_report(reporter_score_res2,rs_threshold = 1.96, y_text_size = 5)
plot2   # Gene_enrichment_Baseline_Tolerance

set.seed(1234)
# Adverse_reaction_occured vs Tolerance_established
# Change of the former (19:26 tolerance established) relative to the latter (10:18 adverse reaction occurred)
reporter_score_res3=reporter_score(KO_abundance[, c(19:26, 10:18)], metadata[c(19:26, 10:18), ], mode="directed", method = "wilcox.test",
                                   p.adjust.method1 = "none", type = "pathway")   #Avoid calculation problems caused by most adjusted p-values being 1.
View(reporter_score_res3$reporter_s)
dim(reporter_score_res3$reporter_s)
plot3 <- plot_report(reporter_score_res3,rs_threshold = 1.96, y_text_size = 5)
plot3   # Gene_enrichment_Intolerance_Tolerance

reporter_score_res <- cbind(reporter_score_res1$reporter_s, reporter_score_res2$reporter_s, reporter_score_res3$reporter_s)
write.csv(reporter_score_res, "results/gene_analysis/reporter_score_results.csv") 

#When we focus on one pathway, e.g. "map02060":
plot_KOs_in_pathway(reporter_score_res1,map_id = "map02060")
plot_KOs_in_pathway(reporter_score_res2,map_id = "map02060")
plot_KOs_in_pathway(reporter_score_res3,map_id = "map02060")

plot_KOs_in_pathway(reporter_score_res1,map_id = "map03440")
plot_KOs_in_pathway(reporter_score_res2,map_id = "map03440")
plot_KOs_in_pathway(reporter_score_res3,map_id = "map03440")

plot_KOs_in_pathway(reporter_score_res1,map_id = "map00280")
plot_KOs_in_pathway(reporter_score_res2,map_id = "map00280")
plot_KOs_in_pathway(reporter_score_res3,map_id = "map00280")

plot_KOs_in_pathway(reporter_score_res1,map_id = "map03050")
plot_KOs_in_pathway(reporter_score_res2,map_id = "map03050")
plot_KOs_in_pathway(reporter_score_res3,map_id = "map03050")

#Or display as a network:
plot_KOs_network(reporter_score_res1,map_id = "map02060", near_pathway = F, main="")
plot_KOs_network(reporter_score_res2,map_id = "map02060", near_pathway = F, main="")
plot_KOs_network(reporter_score_res3,map_id = "map02060", near_pathway = F, main="")

#And we also look at the KOs abundance in a pathway:
plot_KOs_box(reporter_score_res1, map_id = "map02060",only_sig = TRUE)
plot_KOs_box(reporter_score_res2, map_id = "map02060",only_sig = TRUE)
plot_KOs_box(reporter_score_res3, map_id = "map02060",only_sig = TRUE)
## Loading required namespace: ggpubr
## `geom_smooth()` using formula = 'y ~ x'
#Or display as a heatmap:
plot_KOs_heatmap(reporter_score_res1,map_id = "map02060",only_sig = TRUE,heatmap_param = list(cutree_rows=2))
## Loading required namespace: pheatmap


# ############Filter pathways enriched in 2 phases ("Baseline_vs_Intolerance", "Intolerance_vs_Tolerance")
# # that appear enriched (|reporterscore|≥1.96) to create "RS_with_enrich.csv", plot enrichment deviation graph
# #1. Create file

# Step 1: Read
file_path <- "results/gene_analysis/reporter_score_results.csv"
data <- read.csv(file_path, stringsAsFactors = FALSE)
# ---------------------------------------------
# First step: Process first group data (columns 2, 3, 12)
# Intolerance relative to baseline
# ---------------------------------------------
set1 <- data[, c(2, 3, 12)]
colnames(set1) <- c("ID", "Description", "ReporterScore")
set1_filtered <- set1 %>%
  filter(abs(ReporterScore) >= 1.96)
# ---------------------------------------------
# Second step: Process second group data (columns 28, 29, 38)
# Tolerance establishment relative to intolerance occurrence
# ---------------------------------------------
set2 <- data[, c(28, 29, 38)]
colnames(set2) <- c("ID", "Description", "ReporterScore")
set2_filtered <- set2 %>%
  filter(abs(ReporterScore) >= 1.96)
# ---------------------------------------------
# Third step: Find common IDs
# ---------------------------------------------
common_ids <- intersect(set1_filtered$ID, set2_filtered$ID)
# Keep data for common IDs
set1_common <- set1_filtered %>%
  filter(ID %in% common_ids) %>%
  arrange(ID)  # Sort to ensure consistency
set2_common <- set2_filtered %>%
  filter(ID %in% common_ids) %>%
  arrange(ID)  # Sort to ensure consistency
# Check if descriptions for common IDs are consistent
if(!all(set1_common$Description == set2_common$Description)){
  warning("Descriptions for common IDs are inconsistent between the two groups.")
}
# ---------------------------------------------
# Fourth step: Construct result dataframe
# ---------------------------------------------
set1_common <- set1_common %>%
  mutate(ID_Desc = paste(ID, Description, sep = " "))
set2_common <- set2_common %>%
  mutate(ID_Desc = paste(ID, Description, sep = " "))
# Ensure ID_Desc order is consistent between two groups
set1_common <- set1_common %>%
  arrange(ID_Desc)
set2_common <- set2_common %>%
  arrange(ID_Desc)
# Check again if ID_Desc is consistent
if(!all(set1_common$ID_Desc == set2_common$ID_Desc)){
  stop("Common IDs and descriptions are inconsistent in order or content between the two groups.")
}
# Construct result dataframe
result_df <- rbind(
  set1_common$ReporterScore,
  set2_common$ReporterScore
)
colnames(result_df) <- set1_common$ID_Desc
rownames(result_df) <- c("Adverse_reaction_occurred", "Tolerance_established")
result_df <- as.data.frame(result_df)
# ---------------------------------------------
# Fifth step: Save result as CSV file
# ---------------------------------------------
write.csv(result_df, "results/gene_analysis/RS_with_enrich.csv", row.names = TRUE, fileEncoding = "UTF-8")

#Rename the first column to Group

#2. Plotting
pathway_enriched <- read.csv("results/gene_analysis/RS_with_enrich.csv", header =T, check.names = F, sep = ",")
pathway_melt_table <- reshape2::melt(pathway_enriched, id = 1, variable.name = "Pathway", value.name = "Value")
pathway_melt_table$Group <- factor(pathway_melt_table$Group, levels = c("Adverse_reaction_occurred", "Tolerance_established"), ordered = T)
pathway_melt_table$Group
pathway_plot <- ggplot(pathway_melt_table, aes(Value, Pathway, fill=Group, color= Group))+
  geom_col(position = position_dodge2(5,preserve = 'single'), width = 0.5)+
  labs(x="ReporterScore",y=NULL)+
  geom_vline(xintercept = 0)+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size = 1))
pathway_plot    # Enriched gene deviation plot


## Inter-group enriched gene pathway deviation plot
# Prepare KO abundance table and experimental metadata
KO_abundance = read.csv("results/gene_analysis/retinol_KEGG_gene_for_reporter_score.csv", header = T, sep=',', check.names = F, row.names = 1)  #rows are KOs, columns are samples
metadata = read.csv("results/gene_analysis/retinol_group_for_reporter_score.csv", header = T, sep=',', check.names = F, row.names = 1)  #grouping information, first column is sample name, second column is group
load_KOlist()
metadata$Group <- factor(metadata$Group, levels = c("Tolerant_Baseline", "Intolerant_Baseline"), ordered = TRUE)
cat("Comparison: ",levels(factor(metadata$Group)))
#Note: write the compared group (control) at the end, whether the front group (treatment group) is upregulated or downregulated compared to the back group (control group)
#Because the R package has been updated: group as a factor variable, the first level will be set as the control group, you can change the factor level to change your comparison.
#In the code below, control is written at the end, so the reporter score calculation is reversed.
set.seed(1234)
# Change of the former (10:18 intolerant baseline) relative to the latter (1:9 tolerant baseline)
reporter_score_res1=reporter_score(KO_abundance[, c(10:18, 1:9)], metadata[c(10:18 ,1:9), ], mode="directed", method = "wilcox.test",
                                  p.adjust.method1 = "BH", type = "pathway")                             
View(reporter_score_res1$reporter_s)
dim(reporter_score_res1$reporter_s)
plot1 <- plot_report(reporter_score_res1,rs_threshold = 1.96, y_text_size = 5)
plot1   # Gene_enrichment_Baseline_retinol
write.csv(reporter_score_res1$reporter_s, "results/gene_analysis/reporter_score_results_retinol.csv") 

#Manually create RS_with_enrich_retinol.csv

#Plotting
pathway_enriched <- read.csv("results/gene_analysis/RS_with_enrich_retinol.csv", header =T, check.names = F, sep = ",")
pathway_melt_table <- reshape2::melt(pathway_enriched, id = 1, variable.name = "Pathway", value.name = "Value")
pathway_melt_table$Group <- factor(pathway_melt_table$Group, levels = c("Intolerant vs Tolerant"), ordered = T)
pathway_melt_table$Group
pathway_plot <- ggplot(pathway_melt_table, aes(Value, Pathway, fill=Group, color= Group))+
  geom_col(position = position_dodge2(5,preserve = 'single'), width = 0.5)+
  labs(x="ReporterScore",y=NULL)+
  geom_vline(xintercept = 0)+
  geom_vline(xintercept = c(-1.96, 1.96), linetype = "dashed") +  # Add auxiliary lines
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size = 1))
pathway_plot    # Inter-group enriched gene deviation plot