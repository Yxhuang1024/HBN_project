skin_trait <- read.csv("skin_trait.csv", header=T, sep=",", check.names = F) #txt, tsv  #rows as samples, columns as factors
dim(skin_trait)
head(skin_trait)
skin_trait$sample <- factor(skin_trait$sample)
skin_trait$Group <- factor(skin_trait$Group, levels = c("Baseline", "Adverse_reaction_occured", "Tolerance_established"), ordered = T)
skin_trait$Volunteer_ID <- factor(skin_trait$Volunteer_ID)
colnames(skin_trait[,-c(1:3)])[1]

# Group-wise statistical calculations
skin_trait_melt_table <- reshape2::melt(skin_trait, id = 1:3, variable.name = "skin_trait_index", value.name = "Value")
dim(skin_trait_melt_table)
head(skin_trait_melt_table,20)
sink("results/skin_phenotype_analysis/skin_trait_results.csv", split = T)

# Calculate mean and standard deviation
skin_trait_mean_sd_results <- skin_trait_melt_table %>% group_by(skin_trait_index, Group) %>% summarise(mean = mean(Value), sd = sd(Value))
print(skin_trait_mean_sd_results, n=100)

# Statistical significance testing
skin_trait_stat_test <- compare_means(Value ~ Group, group.by = "skin_trait_index", 
                                      paired = T,
                                      data = skin_trait_melt_table)
print(skin_trait_stat_test, n=200)
sink()

# Plotting
comparisons <- list(c("Baseline", "Adverse_reaction_occured"), 
                    c("Adverse_reaction_occured", "Tolerance_established"),
                    c("Baseline", "Tolerance_established"))
skin_trait_violinplot <- ggviolin(skin_trait_melt_table, x = "Group", y = "Value",
                                  color = "Group",fill = "Group",  palette = c("grey","#c97937","#3b6291"), add = "mean_sd",
                                  add.params = list(color = "black", size = 0.2))+
  theme_bw()+
  labs(x = "")+labs(y = "")+
  facet_wrap(~skin_trait_index, ncol = 5, nrow =3, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.8)) +
  stat_compare_means(method = "wilcox.test", paired = TRUE, label = "p.signif", 
                    hide.ns = TRUE, label.y.npc = "top", vjust = 0.5, comparisons = comparisons, p.adjust.method = "BH")

skin_trait_violinplot  # Save as: Three timepoint phenotype significance test


# PCoA plot for all phenotypes
# Environmental factor standardization
data=skin_trait
data[,c(1:3)]
data1=skin_trait[,-c(1:3)]
head(data1)
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

# Construct PCoA table containing PC1 and PC2 axes
species_pcoa <- pcoa_new[, c(2:5)]
head(species_pcoa)
colnames(species_pcoa)[3] <- "PC1"
colnames(species_pcoa)[4] <- "PC2"
dim(species_pcoa)
head(species_pcoa)

# Plotting
species_pcoa$Group <- factor(species_pcoa$Group, levels = c("Baseline", "Adverse_reaction_occured", "Tolerance_established"), ordered = T)
skin_trait_pcoa_plot = ggscatter(species_pcoa, x = "PC1", y = "PC2", 
                                      color = "Group", palette = c("grey","#c97937","#3b6291"), ellipse = TRUE, ellipse.type = "norm",
                                      ellipse.border.remove = F, ellipse.alpha = 0.05, ellipse.level = 0.7, 
                                      size = 1,  legend = "right", mean.point =T, mean.point.size = 4,
                                      xlab = paste0("PCoA1 ( ",pc1,"%"," )",sep="" ),
                                      ylab = paste0("PCoA2 ( ",pc2,"%"," )",sep="" ))+theme_bw()
skin_trait_pcoa_plot  # Save as: PCoA plot of skin phenotype data

# Dissimilarity analysis
sink('results/skin_phenotype_analysis/skin_trait_dissimilarity_test_results.txt', split = T)
set.seed(1234)
pairwise.adonis2(data1~Group, data, nperm = 9999)
sink()    # Save as: skin_trait_dissimilarity_test_results.txt


# PCoA plot for significantly changed phenotypes
# Including phenotypes 1, 3, 4, 8, 11, 12 (need to add 3 to indices below)
data=skin_trait[, c(1:4, 6, 7, 11, 14, 15)] 
data[,c(1:3)]
data1=data[,-c(1:3)]
head(data1)
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

# Construct PCoA table containing PC1 and PC2 axes
species_pcoa <- pcoa_new[, c(2:5)]
head(species_pcoa)
colnames(species_pcoa)[3] <- "PC1"
colnames(species_pcoa)[4] <- "PC2"
dim(species_pcoa)
head(species_pcoa)

# Plotting
species_pcoa$Group <- factor(species_pcoa$Group, levels = c("Baseline", "Adverse_reaction_occured", "Tolerance_established"), ordered = T)
skin_trait_pcoa_plot = ggscatter(species_pcoa, x = "PC1", y = "PC2", 
                                      color = "Group", palette = c("grey","#c97937","#3b6291"), ellipse = TRUE, ellipse.type = "norm",
                                      ellipse.border.remove = F, ellipse.alpha = 0.05, ellipse.level = 0.7, 
                                      size = 1,  legend = "right", mean.point =T, mean.point.size = 4,
                                      xlab = paste0("PCoA1 ( ",pc1,"%"," )",sep="" ),
                                      ylab = paste0("PCoA2 ( ",pc2,"%"," )",sep="" ))+theme_bw()
skin_trait_pcoa_plot  # Save as: PCoA plot of significantly changed skin phenotype data

# Dissimilarity analysis
sink('results/skin_phenotype_analysis/significant_skin_trait_dissimilarity_test_results.txt', split = T)
set.seed(1234)
pairwise.adonis2(data1~Group, data, nperm = 9999)
sink()    # Save as: significant_skin_trait_dissimilarity_test_results.txt



####### Phenotype differences between 9 retinol-tolerant and 9 intolerant volunteers #######
skin_trait <- read.csv("skin_trait_retinol_yes_no.csv", header=T, sep=",", check.names = F) #txt, tsv  #rows as samples, columns as factors
dim(skin_trait)
head(skin_trait)
skin_trait$sample <- factor(skin_trait$sample)
skin_trait$Group <- factor(skin_trait$Group, levels = c("Intolerant_Baseline", "Tolerant_Baseline"), ordered = T)
colnames(skin_trait[,-c(1:2)])[1]

# Group-wise statistical calculations
skin_trait_melt_table <- reshape2::melt(skin_trait, id = 1:2, variable.name = "skin_trait_index", value.name = "Value")
dim(skin_trait_melt_table)
head(skin_trait_melt_table,20)
sink("results/skin_phenotype_analysis/skin_trait_retinol_yes_no_results.csv", split = T)

# Calculate mean and standard deviation
skin_trait_mean_sd_results <- skin_trait_melt_table %>% group_by(skin_trait_index, Group) %>% summarise(mean = mean(Value), sd = sd(Value))
print(skin_trait_mean_sd_results, n=100)

# Statistical significance testing
skin_trait_stat_test <- compare_means(Value ~ Group, group.by = "skin_trait_index", 
                                      #paired = T,
                                      data = skin_trait_melt_table)
print(skin_trait_stat_test, n=200)
sink()

# Plotting
comparisons <- list(c("Intolerant_Baseline", "Tolerant_Baseline"))
skin_trait_violinplot <- ggviolin(skin_trait_melt_table, x = "Group", y = "Value",
                                  color = "Group",fill = "Group",  palette = c("#c97937","#3b6291"), add = "mean_sd",
                                  add.params = list(color = "black", size = 0.2))+
  theme_bw()+
  labs(x = "")+labs(y = "")+
  facet_wrap(~skin_trait_index, ncol = 5, nrow =3, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.8)) +
  stat_compare_means(method = "wilcox.test", paired = FALSE, label = "p.signif", 
                    hide.ns = FALSE, label.y.npc = "top", comparisons = comparisons, p.adjust.method = "BH")
skin_trait_violinplot  # Save as: Inter-group baseline skin phenotype comparison

# PCoA plot
# Environmental factor standardization
data=skin_trait
data[,c(1:2)]
data1=skin_trait[,-c(1:2)]
head(data1)
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

# Construct PCoA table containing PC1 and PC2 axes
species_pcoa <- pcoa_new[, c(2:5)]
head(species_pcoa)
colnames(species_pcoa)[3] <- "PC1"
colnames(species_pcoa)[4] <- "PC2"
dim(species_pcoa)
head(species_pcoa)

# Plotting
species_pcoa$Group <- factor(species_pcoa$Group, levels = c("Intolerant_Baseline", "Tolerant_Baseline"), ordered = T)
skin_trait_pcoa_plot = ggscatter(species_pcoa, x = "PC1", y = "PC2", 
                                 color = "Group", palette = c("#c97937","#3b6291"), ellipse = TRUE, ellipse.type = "norm",
                                 ellipse.border.remove = F, ellipse.alpha = 0.05, ellipse.level = 0.7, 
                                 size = 1,  legend = "right", mean.point =T, mean.point.size = 4,
                                 xlab = paste0("PCoA1 ( ",pc1,"%"," )",sep="" ),
                                 ylab = paste0("PCoA2 ( ",pc2,"%"," )",sep="" ))+theme_bw()
skin_trait_pcoa_plot  # Save as: Inter-group baseline skin phenotype comparison PCoA plot

# Dissimilarity analysis
sink('results/skin_phenotype_analysis/skin_trait_retinol_yes_no_dissimilarity_test_results.txt', split = T)
set.seed(1234)
pairwise.adonis2(data1~Group, data, nperm = 9999)
sink()