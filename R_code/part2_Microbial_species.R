##############     Calculate species/OTU diversity     ##############
species <- read.csv("mOTUs3_species.txt", header =T, check.names = F, sep = "\t")
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
write.csv(shannon_table, "results/species_analysis/species_shannon_index_table_all_species.csv")
dim(shannon_table)
head(shannon_table)
sink("results/species_analysis/species_shannon_results.csv", split = T)
#calculate mean_sd
shannon_results <- shannon_table %>% group_by(Group) %>% summarise(mean = mean(shannon), sd = sd(shannon))
print(shannon_results)
#calculate significance
shannon_stat_test <- compare_means(shannon ~ Group, 
                                   #paired = T, 
                                   data = shannon_table)
print(shannon_stat_test)
sink()

#richness
richness <-as.matrix(specnumber(species[,-c(1:3)]))
richness
richness_table <- cbind.data.frame(species[,1:3], richness)
richness_table
write.csv(richness_table, "results/species_analysis/species_richness_table_all_species.csv") 
dim(richness_table)
head(richness_table)
sink("results/species_analysis/species_richness_results.csv", split = T)
#calculate mean_sd
richness_results <- richness_table %>% group_by(Group) %>% summarise(mean = mean(richness), sd = sd(richness))
print(richness_results)
#calculate significance
richness_stat_test <- compare_means(richness ~ Group, 
                                    #paired = T, 
                                    data = richness_table)
print(richness_stat_test)
sink()

#evenness
evenness <- as.matrix(diversity(species[,-c(1:3)],index="shannon")/log(specnumber(species[,-c(1:3)])))
evenness
evenness_table <- cbind.data.frame(species[,1:3], evenness)
evenness_table
write.csv(evenness_table, "results/species_analysis/species_evenness_table_all_species.csv") 
dim(evenness_table)
head(evenness_table)
sink("results/species_analysis/species_evenness_results.csv", split = T)
#calculate mean_sd
evenness_results <- evenness_table %>% group_by(Group) %>% summarise(mean = mean(evenness), sd = sd(evenness))
print(evenness_results)
#calculate significance
evenness_stat_test <- compare_means(evenness ~ Group, 
                                    #paired = T, 
                                    data = evenness_table)
print(evenness_stat_test)
sink()

#plotting
comparisons <- list(c("Baseline", "Adverse_reaction_occured"), 
                    c("Adverse_reaction_occured", "Tolerance_established"), 
                    c("Baseline", "Tolerance_established"))
# Shannon violin plot
shannon_violinplot <- ggviolin(shannon_table, x = "Group", y = "shannon",
                               color = "Group", fill = "Group", palette = c("grey", "#c97937", "#3b6291"), add = "mean_sd",
                               add.params = list(color = "white", size = 0.5)) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif", paired = FALSE) +
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
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif", paired = FALSE) +
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
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif", paired = FALSE) +
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
alpha_diversity_violinplot  # Save as: Species α diversity

#pcoa plot
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
species_pcoa_plot   # Save as: Species PCoA plot
######## adonis dissimilarity analysis
sink('results/species_analysis/species_dissimilarity_test_results.txt', split = T)
set.seed(1234)
pairwise.adonis2(data1~Group, data, nperm = 9999)
sink()

#Statistical analysis of TOP20 abundance and biological classification species information
species <- read.csv("mOTUs3_species.txt", header = TRUE, check.names = FALSE, sep = "\t")
taxonomy_info <- colnames(species)[4:ncol(species)]
taxonomy_levels <- strsplit(taxonomy_info, split = "\\|")
taxonomy_df <- do.call(rbind, taxonomy_levels)
colnames(taxonomy_df) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxonomy_count <- apply(taxonomy_df, 2, function(x) length(unique(x)))
print(taxonomy_count)
high_frequency_species <- sort(table(taxonomy_df[, "Species"]), decreasing = TRUE)
print(head(high_frequency_species, 10))
abundance_data <- species[, 4:ncol(species)]
average_abundance <- colMeans(abundance_data, na.rm = TRUE)
top_20_species <- sort(average_abundance, decreasing = TRUE)[1:20]
top_20_species_names <- names(top_20_species)
top_20_species_df <- data.frame(Species = top_20_species_names, Avg_Abundance = top_20_species)
print(top_20_species_df)
write.csv(taxonomy_count, "taxonomy_count_R.csv", row.names = FALSE)
write.csv(head(high_frequency_species, 10), "high_frequency_species_R.csv", row.names = FALSE)
write.csv(top_20_species_df, "top_20_species_R.csv", row.names = FALSE)

#Top20 species composition plot (each sample)
species_top20 <- read.csv("mOTUs3_species_top20.csv", header=T, sep=",", check.names = F)
head(species_top20)
melt_data <- reshape2::melt(species_top20, id = 1:3, variable.name = "Species", value.name = "Relative_abundance")
dim(melt_data)
head(melt_data, 20)

col <- c("#5AD1E7","#A1BA28","#C5C6DE","#D75CEA","#E0CC9D",
                  "#F9E7E9","firebrick","#43D0AC","#6624B4","#6EB7AF",
                  "#F762B1","#E1FDA4","#84BBED","#DCC063","#F3F862",
                  "#89806A","#EEC2B6","#EE819D","#9A7CDF","#D1F9FE","darkgrey")  #top20----22 colors
                  
melt_data$Group <- factor(melt_data$Group, levels = c("Baseline", "Adverse_reaction_occured", "Tolerance_established"), ordered = TRUE)
melt_data$Group
melt_data$Volunteer_ID <- factor(melt_data$Volunteer_ID, ordered = TRUE)
melt_data$Volunteer_ID 

mOTUs3_top20_species_plot <- ggbarplot(melt_data, x = "Volunteer_ID", y = "Relative_abundance", 
                                       xlab = paste0("Volunteer ID"),
                                       color = "Species", fill = "Species", palette = col, 
                                       legend = "top", scales = "free_x", 
                                       ylab = paste0("Relative abundance (%)")) +
  facet_wrap(~ Group, ncol = 1) +  # Arrange plots vertically
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 9, face = "italic", family = "sans")) +
  guides(fill = guide_legend(ncol = 1)) +  # Display legend in one column
  theme(legend.position = "right") +  # Legend on the right
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +  # Add plot border
  ggtitle("Top20 species") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

mOTUs3_top20_species_plot  # Save as: TOP20 species

#Top20 species composition plot (by group)
species_top20 <- read.csv("mOTUs3_species_top20.csv", header = TRUE, sep = ",", check.names = FALSE)
melt_data <- reshape2::melt(species_top20, id = 1:3, variable.name = "Species", value.name = "Relative_abundance")
col <- c("#5AD1E7", "#A1BA28", "#C5C6DE", "#D75CEA", "#E0CC9D", 
         "#F9E7E9", "firebrick", "#43D0AC", "#6624B4", "#6EB7AF", 
         "#F762B1", "#E1FDA4", "#84BBED", "#DCC063", "#F3F862", 
         "#89806A", "#EEC2B6", "#EE819D", "#9A7CDF", "#D1F9FE", "darkgrey")
melt_data$Group <- factor(melt_data$Group, levels = c("Baseline", "Adverse_reaction_occured", "Tolerance_established"), ordered = TRUE)
melt_data$Volunteer_ID <- factor(melt_data$Volunteer_ID, ordered = TRUE)
# Calculate average relative abundance for each species within each group
group_avg_data <- melt_data %>%
  group_by(Group, Species) %>%
  summarise(Average_Relative_abundance = mean(Relative_abundance, na.rm = TRUE)) %>%
  ungroup()
# View first few rows of the dataframe
head(group_avg_data)
# Create stacked plot
group_avg_plot <- ggplot(group_avg_data, aes(x = Group, y = Average_Relative_abundance, fill = Species)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = col) +  # Use your defined colors
  labs(x = "Group", y = "Average Relative Abundance (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 9, face = "italic", family = "sans")) +
  guides(fill = guide_legend(ncol = 1)) +  # Display legend in one column
  theme(legend.position = "right") +  # Legend on the right
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +  # Add plot border
  ggtitle("Average Top 20 Species Abundance by Group") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
group_avg_plot  # Save as: TOP20 species (by group)

#Significance testing for species
mOTUs3_all_species <- read.csv("mOTUs3_species.txt", header = TRUE, check.names = FALSE, sep = "\t")
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
write.csv(species_abundance_result, "results/species_analysis/species_abundance_mean_sd_Group.csv")
data1 <- species_melt_table
species_stat_test1 <- compare_means(Value ~ Group, group.by = c("Species"), paired = FALSE, data = data1)
print(species_stat_test1, n = 200)
write.csv(species_stat_test1, "results/species_analysis/species_abundance_stat_test_Group.csv")

#Plots for significantly changed species
#Removed all-zero data: mOTUs3_species_with_signif_remove_0.txt
# skin_trait <- read.csv("mOTUs3_species_with_signif_remove_0.txt", header = TRUE, check.names = FALSE, sep = "\t")
skin_trait <- read.csv("mOTUs3_species_with_signif.txt", header = TRUE, check.names = FALSE, sep = "\t")
dim(skin_trait)
head(skin_trait)
skin_trait$sample <- factor(skin_trait$sample)
skin_trait$Group <- factor(skin_trait$Group, levels = c("Baseline", "Adverse_reaction_occured", "Tolerance_established"), ordered = T)
skin_trait$Volunteer_ID <- factor(skin_trait$Volunteer_ID)
colnames(skin_trait[,-c(1:3)])[1]
#Group-wise statistical calculations
skin_trait_melt_table <- reshape2::melt(skin_trait, id = 1:3, variable.name = "skin_trait_index", value.name = "Value")
dim(skin_trait_melt_table)
head(skin_trait_melt_table,20)
#plotting
comparisons <- list(c("Baseline", "Adverse_reaction_occured"), 
                    c("Adverse_reaction_occured", "Tolerance_established"),
                    c("Baseline", "Tolerance_established"))
# Remove square brackets and their contents
skin_trait_melt_table$skin_trait_index <- gsub("\\[.*\\]", "", skin_trait_melt_table$skin_trait_index)
skin_trait_violinplot <- ggboxplot(skin_trait_melt_table, x = "Group", y = "Value",
                                  color = "Group",fill = "Group",  palette = c("grey","#c97937","#3b6291"), add = "mean_sd",
                                  add.params = list(color = "black", size = 0.2))+
  theme_bw()+
  labs(x = "")+labs(y = "")+
  facet_wrap(~skin_trait_index, ncol = 6, nrow =3, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.8)) +
  theme(strip.text = element_text(size = 8)) +  # Adjust subtitle font size
  stat_compare_means(method = "wilcox.test", paired = FALSE, label = "p.signif", 
                    hide.ns = TRUE, label.y.npc = "bottom", vjust = 0.5,  comparisons = comparisons, p.adjust.method = "BH")

skin_trait_violinplot  # Save as: Three timepoint species significance testing


# Plot significance testing for these microorganisms at three timepoints in tolerance group (17 significantly changed microorganisms mentioned above)
# Modified complete code
skin_trait <- read.csv("results/species_analysis/mOTUs3_species_tolerance_group_17species.csv", header = TRUE, check.names = FALSE, sep = ",")
# Check original timepoints
print("Original unique values in Time column:")
print(unique(skin_trait$Time))
# Keep only needed timepoints
skin_trait <- skin_trait[skin_trait$Time %in% c("Day0", "Day7", "Day21"), ]
# Check timepoints after filtering
print("Unique values in Time column after filtering:")
print(unique(skin_trait$Time))
# Data preprocessing
skin_trait$sample <- factor(skin_trait$sample)
skin_trait$Time <- factor(skin_trait$Time, levels = c("Day0", "Day7", "Day21"), ordered = T)
# Reshape data
skin_trait_melt_table <- reshape2::melt(skin_trait, id = 1:2, variable.name = "skin_trait_index", value.name = "Value")
# Set comparison groups
comparisons <- list(c("Day0", "Day7"), 
                    c("Day7", "Day21"),
                    c("Day0", "Day21"))
# Remove square brackets and their contents
skin_trait_melt_table$skin_trait_index <- gsub("\\[.*\\]", "", skin_trait_melt_table$skin_trait_index)
# Create boxplot
skin_trait_violinplot <- ggboxplot(skin_trait_melt_table, x = "Time", y = "Value",
                                  color = "Time", fill = "Time", 
                                  palette = c("#2E8B57", "#FF6347", "#4682B4"), 
                                  add = "mean_sd",
                                  add.params = list(color = "black", size = 0.2)) +
  theme_bw() +
  labs(x = "", y = "") +
  facet_wrap(~skin_trait_index, ncol = 5, nrow = 3, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.8)) +
  theme(strip.text = element_text(size = 8)) +
  stat_compare_means(method = "wilcox.test", 
                    paired = FALSE,  # Changed to unpaired test
                    label = "p.signif", 
                    hide.ns = FALSE,  # Show all results
                    comparisons = comparisons, 
                    p.adjust.method = "BH",
                    label.y.npc = 0,
                    step.increase = 0.1,  # Adjust label spacing
                    tip.length = 0.01,  # Adjust connecting line length
                    vjust = 0.25,
                    size = 3)    

skin_trait_violinplot  # Save as: 17 significantly changed species at four timepoints in tolerance group species significance testing
                       # sup fig 2



#Plots for species that significantly changed in both phases
#Including 1. tolerance establishment relative to adverse reaction occurrence;
#         2. adverse reaction occurrence relative to baseline
skin_trait <- read.csv("mOTUs3_species_with_signif.txt", header = TRUE, check.names = FALSE, sep = "\t")
dim(skin_trait)
head(skin_trait)
skin_trait$sample <- factor(skin_trait$sample)
skin_trait$Group <- factor(skin_trait$Group, levels = c("Baseline", "Adverse_reaction_occured", "Tolerance_established"), ordered = T)
skin_trait$Volunteer_ID <- factor(skin_trait$Volunteer_ID)
colnames(skin_trait[,-c(1:3)])[1]
#Group-wise statistical calculations
skin_trait_melt_table <- reshape2::melt(skin_trait, id = 1:3, variable.name = "skin_trait_index", value.name = "Value")
dim(skin_trait_melt_table)
head(skin_trait_melt_table,20)
#plotting
comparisons <- list(c("Baseline", "Adverse_reaction_occured"), 
                    c("Adverse_reaction_occured", "Tolerance_established"),
                    c("Baseline", "Tolerance_established"))
# Remove square brackets and their contents
skin_trait_melt_table$skin_trait_index <- gsub("\\[.*\\]", "", skin_trait_melt_table$skin_trait_index)
skin_trait_violinplot <- ggboxplot(skin_trait_melt_table, x = "Group", y = "Value",
                                  color = "Group",fill = "Group",  palette = c("grey","#c97937","#3b6291"), add = "mean_sd",
                                  add.params = list(color = "black", size = 0.2))+
  theme_bw()+
  labs(x = "")+labs(y = "")+
  facet_wrap(~skin_trait_index, ncol = 6, nrow =3, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.8)) +
  theme(strip.text = element_text(size = 8)) +  # Adjust subtitle font size
  stat_compare_means(method = "wilcox.test", paired = FALSE, label = "p.signif", 
                    hide.ns = TRUE, label.y.npc = "top", vjust = 0.5,  comparisons = comparisons, p.adjust.method = "BH")
skin_trait_violinplot  # Save as: Three timepoint species significance testing (changed in both phases)

#Save data for species that significantly changed in both phases


#PCoA plot for significantly changed species
species <- read.csv("mOTUs3_species_with_signif.txt", header = TRUE, check.names = FALSE, sep = "\t")
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
species_pcoa_plot   # Save as: PCoA plot of significantly changed species
######## adonis dissimilarity analysis
sink('results/species_analysis/significant_species_dissimilarity_test_results.txt', split = T)
set.seed(1234)
pairwise.adonis2(data1~Group, data, nperm = 9999)
sink()



####### Significantly changed microbial species differences between 9 retinol-tolerant and 9 intolerant volunteers #######
## Create species file
# See get_compare_data.py

## α diversity comparison (can only plot violin plots, line plots cannot be drawn due to unpaired data (ID))
species <- read.csv("results/species_analysis/mOTUs3_species_retinol_yes_no_new.csv", header=T, sep=",", check.names = F) 
dim(species)
head(species)
# 1. Ensure relevant columns are factors
species$sample <- factor(species$sample)
species$Group <- factor(species$Group, levels = c("Intolerant_Baseline", "Tolerant_Baseline")) # Adjust according to actual group names
species$Volunteer_ID <- factor(species$Volunteer_ID)
# 2. Calculate Shannon index
shannon <- diversity(species[,-c(1:2)], index = "shannon")
shannon_table <- cbind.data.frame(species[,1:2], shannon)
write.csv(shannon_table, "results/species_analysis/species_shannon_index_table_retinol_yes_no.csv", row.names = FALSE)
dim(shannon_table)
head(shannon_table)
# Statistics and significance testing
sink("results/species_analysis/species_shannon_results_retinol_yes_no.csv", split = TRUE)
# Calculate mean and sd
shannon_results <- shannon_table %>% group_by(Group) %>% summarise(mean = mean(shannon), sd = sd(shannon))
print(shannon_results)
# Calculate significance (Wilcoxon test, unpaired)
shannon_stat_test <- compare_means(shannon ~ Group, 
                                   method = "wilcox.test", 
                                   data = shannon_table)
print(shannon_stat_test)
sink()
# 3. Calculate Richness
richness <- specnumber(species[,-c(1:2)])
richness_table <- cbind.data.frame(species[,1:2], richness)
write.csv(richness_table, "results/species_analysis/species_richness_table_retinol_yes_no.csv", row.names = FALSE)
dim(richness_table)
head(richness_table)
# Statistics and significance testing
sink("results/species_analysis/species_richness_results_retinol_yes_no.csv", split = TRUE)
# Calculate mean and sd
richness_results <- richness_table %>% group_by(Group) %>% summarise(mean = mean(richness), sd = sd(richness))
print(richness_results)
# Calculate significance (Wilcoxon test, unpaired)
richness_stat_test <- compare_means(richness ~ Group, 
                                    method = "wilcox.test", 
                                    data = richness_table)
print(richness_stat_test)
sink()
# 4. Calculate Evenness
evenness <- diversity(species[,-c(1:2)], index = "shannon") / log(specnumber(species[,-c(1:2)]))
evenness_table <- cbind.data.frame(species[,1:2], evenness)
write.csv(evenness_table, "results/species_analysis/species_evenness_table_retinol_yes_no.csv", row.names = FALSE)
dim(evenness_table)
head(evenness_table)
# Statistics and significance testing
sink("results/species_analysis/species_evenness_results_retinol_yes_no.csv", split = TRUE)
# Calculate mean and sd
evenness_results <- evenness_table %>% group_by(Group) %>% summarise(mean = mean(evenness), sd = sd(evenness))
print(evenness_results)
# Calculate significance (Wilcoxon test, unpaired)
evenness_stat_test <- compare_means(evenness ~ Group, 
                                    method = "wilcox.test", 
                                    data = evenness_table)
print(evenness_stat_test)
sink()
# 5. Plot α diversity
# Define comparison groups
comparisons <- list(c("Intolerant_Baseline", "Tolerant_Baseline"))
# 5.1 Shannon violin plot
shannon_violinplot <- ggviolin(shannon_table, x = "Group", y = "shannon",
                               color = "Group", fill = "Group", palette = c("#c97937", "#3b6291"), add = "mean_sd",
                               add.params = list(color = "white", size = 0.5)) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif", paired = FALSE) +
  labs(y = "Species Shannon index", x = "") +
  theme_bw()
# 5.3 Richness violin plot
richness_violinplot <- ggviolin(richness_table, x = "Group", y = "richness",
                                color = "Group", fill = "Group", palette = c("#c97937", "#3b6291"), add = "mean_sd",
                                add.params = list(color = "white", size = 0.5)) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif", paired = FALSE) +
  labs(y = "Species Richness", x = "") +
  theme_bw()
# 5.5 Evenness violin plot
evenness_violinplot <- ggviolin(evenness_table, x = "Group", y = "evenness",
                                color = "Group", fill = "Group", palette = c("#c97937", "#3b6291"), add = "mean_sd",
                                add.params = list(color = "white", size = 0.5)) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif", paired = FALSE) +
  labs(y = "Species Evenness", x = "") +
  theme_bw()
# 5.7 Combine α diversity plots
alpha_diversity_violinplot <- ggarrange(
  ggarrange(shannon_violinplot, richness_violinplot, evenness_violinplot, 
            nrow = 3, ncol = 1, common.legend = TRUE, legend = "right"),
  nrow = 1, ncol = 2
)   
alpha_diversity_violinplot  # Save as: results/species_analysis/alpha_diversity_retinol_yes_no


## Baseline TOP20 species comparison between two groups
## 1. Read data
species_data <- read.csv("results/species_analysis/mOTUs3_species_retinol_yes_no_new.csv", 
                           header = TRUE, sep = ",", check.names = FALSE)

## 2. Calculate TOP20 species for each group
compute_top20 <- function(df, group_name) {
  df_group <- df %>% filter(Group == group_name)
  df_group_long <- df_group %>%
    select(-sample, -Group) %>%
    pivot_longer(cols = everything(), names_to = "Species", values_to = "Relative_abundance")
  
  top20_species <- df_group_long %>%
    group_by(Species) %>%
    summarise(mean_abundance = mean(Relative_abundance, na.rm = TRUE)) %>%
    arrange(desc(mean_abundance)) %>%
    slice_head(n = 20) %>%
    pull(Species)
  
  return(top20_species)
}

top20_intolerant_list <- compute_top20(species_data, "Intolerant_Baseline")
top20_tolerant_list   <- compute_top20(species_data, "Tolerant_Baseline")

## 3. Aggregate TOP20 data, classify non-TOP20 as Others, and normalize to calculate percentages
aggregate_top20 <- function(df, group_name, top_list) {
  df_group <- df %>% filter(Group == group_name)
  df_group_long <- df_group %>%
    select(-sample, -Group) %>%
    pivot_longer(cols = everything(), names_to = "Species", values_to = "Relative_abundance")
  
  # Classify non-TOP20 as "Others"
  df_group_long <- df_group_long %>%
    mutate(Species = ifelse(Species %in% top_list, Species, "Others"))
  
  # Aggregate abundance for same species and calculate normalized percentages
  df_group_summary <- df_group_long %>%
    group_by(Species) %>%
    summarise(Relative_abundance = sum(Relative_abundance, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(Group = group_name) %>%
    group_by(Group) %>%
    mutate(Total = sum(Relative_abundance),
           Relative_abundance_percent = (Relative_abundance / Total) * 100) %>%
    ungroup()
  
  return(df_group_summary)
}

top20_intolerant_df <- aggregate_top20(species_data, "Intolerant_Baseline", top20_intolerant_list)
top20_tolerant_df   <- aggregate_top20(species_data, "Tolerant_Baseline", top20_tolerant_list)

## 4. Save intermediate results
output_dir <- "results/species_analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
write.csv(top20_intolerant_df, file = file.path(output_dir, "Intolerant_Top20_species.csv"), row.names = FALSE)
write.csv(top20_tolerant_df, file = file.path(output_dir, "Tolerant_Top20_species.csv"), row.names = FALSE)

## 5. Combine data from both groups, construct unified color system
combined_df <- bind_rows(top20_intolerant_df, top20_tolerant_df)

# Define factor order: sort non-Others species first, then put Others last
species_order <- combined_df %>% 
  filter(Species != "Others") %>% 
  group_by(Species) %>% 
  summarise(total = sum(Relative_abundance, na.rm = TRUE)) %>%
  arrange(desc(total)) %>%
  pull(Species)
common_levels <- c(species_order, "Others")
combined_df$Species <- factor(combined_df$Species, levels = common_levels)

# Define color palette
user_colors <- c("#5AD1E7", "#A1BA28", "#C5C6DE", "#D75CEA", "#E0CC9D", 
                 "#F9E7E9", "firebrick", "#43D0AC", "#6624B4", "#6EB7AF", 
                 "#F762B1", "#E1FDA4", "#84BBED", "#DCC063", "#F3F862", 
                 "#89806A", "#EEC2B6", "#EE819D", "#9A7CDF", "#D1F9FE", "darkgrey")
non_grey_colors <- user_colors[user_colors != "darkgrey"]

n_species <- length(species_order)
if(n_species > length(non_grey_colors)){
  non_grey_colors <- colorRampPalette(non_grey_colors)(n_species)
}
common_palette <- setNames(non_grey_colors[1:n_species], species_order)
common_palette["Others"] <- "darkgrey"

## 6. Create combined stacked bar plot (using normalized percentage data)
combined_plot <- ggplot(combined_df, aes(x = Group, y = Relative_abundance_percent, fill = Species)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5) +
  scale_fill_manual(values = common_palette, guide = guide_legend(ncol = 1)) +
  labs(x = "Group", y = "Relative Abundance (%)", fill = "Species") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 9, face = "italic"),
    legend.key.size = unit(0.4, "cm"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  ) +
  ggtitle("Relative Abundance of Top 20 Species by Group")

## 7. Display plot
print(combined_plot)




## Differential species comparison
skin_trait <- read.csv("results/species_analysis/mOTUs3_species_retinol_yes_no_new.csv", header=T, sep=",", check.names = F) 
dim(skin_trait)
head(skin_trait)
skin_trait$sample <- factor(skin_trait$sample)
skin_trait$Group <- factor(skin_trait$Group, levels = c("Intolerant_Baseline", "Tolerant_Baseline"), ordered = T)
colnames(skin_trait[,-c(1:2)])[1]
#Group-wise statistical calculations
skin_trait_melt_table <- reshape2::melt(skin_trait, id = 1:2, variable.name = "skin_trait_index", value.name = "Value")
dim(skin_trait_melt_table)
head(skin_trait_melt_table,18)
#calculate mean_sd
skin_trait_mean_sd_results <- skin_trait_melt_table %>% group_by(skin_trait_index, Group) %>% summarise(mean = mean(Value), sd = sd(Value))
write.csv(skin_trait_mean_sd_results, "results/species_analysis/species_retinol_yes_no_mean_sd_results.csv", row.names = FALSE)
#calculate significance
skin_trait_stat_test <- compare_means(Value ~ Group, group.by = "skin_trait_index", 
                                      #paired = T,
                                      data = skin_trait_melt_table,
                                      p.adjust.method = "BH")
write.csv(skin_trait_stat_test, "results/species_analysis/species_retinol_yes_no_stat_test.csv", row.names = FALSE)

#plotting                                  
skin_trait <- read.csv("results/species_analysis/mOTUs3_species_retinol_yes_no_new.csv", header = TRUE, sep = ",", check.names = FALSE)
# Data preprocessing
skin_trait$sample <- factor(skin_trait$sample)
skin_trait$Group <- factor(skin_trait$Group, levels = c("Intolerant_Baseline", "Tolerant_Baseline"), ordered = TRUE)
# Determine columns to plot, here you need to manually find the positions of these columns according to the above file and rename these columns
selected_cols <- colnames(skin_trait)[c(176, 256, 338, 361, 445, 1063, 1139, 1331, 1336, 1343, 1414)]
# Data transformation (melt)
skin_trait_melt_table <- melt(skin_trait, id = c("sample", "Group"), variable.name = "skin_trait_index", value.name = "Value")
# Filter out specified columns
skin_trait_melt_table_selected <- skin_trait_melt_table %>%
  filter(skin_trait_index %in% selected_cols)
# plotting
comparisons <- list(c("Intolerant_Baseline", "Tolerant_Baseline"))
skin_trait_violinplot <- ggboxplot(skin_trait_melt_table_selected, x = "Group", y = "Value",
                                  color = "Group", fill = "Group", palette = c("#c97937", "#3b6291"),
                                  add = "mean_sd",
                                  add.params = list(color = "black", size = 0.2)) +
  theme_bw() +
  labs(x = "", y = "") +
  facet_wrap(~skin_trait_index, ncol = 4, nrow = 3, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.8)) +
  stat_compare_means(method = "wilcox.test", paired = FALSE, label = "p.signif", 
                     hide.ns = FALSE, label.y.npc = "top", vjust = 0.5, 
                     comparisons = comparisons, p.adjust.method = "BH")
skin_trait_violinplot # Inter-group baseline species comparison

#pcoa plot
#Environmental factor standardization
skin_trait <- read.csv("results/species_analysis/mOTUs3_species_retinol_yes_no_new.csv", header=T, sep=",", check.names = F) 
row_sums <- rowSums(skin_trait[, -c(1, 2)] )
skin_trait_filtered <- skin_trait[row_sums != 0, ]  # Filter all-zero rows
skin_trait <- na.omit(skin_trait_filtered)  # Filter missing values

data=skin_trait
data1=skin_trait[,-c(1:2)]
head(data1)
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
skin_trait_pcoa_plot  # Save as: Inter-group baseline species comparison PCoA plot
#Dissimilarity analysis
sink('results/species_analysis/species_retinol_yes_no_dissimilarity_test_results.txt', split = T)
set.seed(1234)
pairwise.adonis2(data1~Group, data, nperm = 9999)
sink()