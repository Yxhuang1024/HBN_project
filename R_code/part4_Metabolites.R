###############     metabolites analysis    ###########################
#Replace metabolite abundance 0.0000099 with 0
species <- read.csv("metabolites.csv", header =T, check.names = F, sep = ",")
species[species == 0.0000099] <- 0
write.csv(species, "metabolites.csv", row.names = FALSE)

#Read data
species <- read.csv("metabolites.csv", header =T, check.names = F, sep = ",")
species$sample <- factor(species$sample)
species$Group <- factor(species$Group, levels = c("Baseline", "Adverse_reaction_occurred", "Tolerance_established"))
species$Volunteer_ID <- factor(species$Volunteer_ID)
colnames(species[,-c(1:3)])[1]

#shannon
shannon <-as.matrix(diversity(species[,-c(1:3)],index="shannon"))
shannon
shannon_table <- cbind.data.frame(species[,1:3], shannon)
shannon_table
write.csv(shannon_table, "results/metabolite_analysis/metabolites_shannon_index_table.csv")
dim(shannon_table)
head(shannon_table)
sink("results/metabolite_analysis/metabolites_shannon_results.csv", split = T)
#calculate mean_sd
shannon_results <- shannon_table %>% group_by(Group) %>% summarise(mean = mean(shannon), sd = sd(shannon))
print(shannon_results)
#calculate significance
shannon_stat_test <- compare_means(shannon ~ Group, 
                                   paired = TRUE, 
                                   data = shannon_table,
                                   p.adjust.method1 = "BH")
print(shannon_stat_test)
sink()

#richness
richness <-as.matrix(specnumber(species[,-c(1:3)]))
richness
richness_table <- cbind.data.frame(species[,1:3], richness)
richness_table
write.csv(richness_table, "results/metabolite_analysis/metabolites_richness_table.csv") 
dim(richness_table)
head(richness_table)
sink("results/metabolite_analysis/metabolites_richness_results.csv", split = T)
#calculate mean_sd
richness_results <- richness_table %>% group_by(Group) %>% summarise(mean = mean(richness), sd = sd(richness))
print(richness_results)
#calculate significance
richness_stat_test <- compare_means(richness ~ Group, 
                                    paired = TRUE, 
                                    data = richness_table,
                                    p.adjust.method1 = "BH")
print(richness_stat_test)
sink()

#evenness
evenness <- as.matrix(diversity(species[,-c(1:3)],index="shannon")/log(specnumber(species[,-c(1:3)])))
evenness
evenness_table <- cbind.data.frame(species[,1:3], evenness)
evenness_table
write.csv(evenness_table, "results/metabolite_analysis/metabolites_evenness_table.csv") 
dim(evenness_table)
head(evenness_table)
sink("results/metabolite_analysis/metabolites_evenness_results.csv", split = T)
#calculate mean_sd
evenness_results <- evenness_table %>% group_by(Group) %>% summarise(mean = mean(evenness), sd = sd(evenness))
print(evenness_results)
#calculate significance
evenness_stat_test <- compare_means(evenness ~ Group, 
                                    paired = TRUE,  
                                    data = evenness_table,
                                    p.adjust.method1 = "BH")
print(evenness_stat_test)
sink()

#plotting
comparisons <- list(c("Baseline", "Adverse_reaction_occurred"), 
                    c("Adverse_reaction_occurred", "Tolerance_established"), 
                    c("Baseline", "Tolerance_established"))
# Shannon violin plot
shannon_violinplot <- ggviolin(shannon_table, x = "Group", y = "shannon",
                               color = "Group", fill = "Group", palette = c("grey", "#c97937", "#3b6291"), add = "mean_sd",
                               add.params = list(color = "white", size = 0.5)) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif", paired = TRUE, p.adjust.method1 = "BH") +
  labs(y = "metabolites Shannon index", x = "") +
  theme_bw()
shannon_violinplot
# Shannon line plot
shannon_lineplot <- ggline(shannon_table, x = "Group", y = "shannon", group = "Volunteer_ID",
                           size = 0.3, color = "Volunteer_ID", palette = "lancet") +
  labs(x = NULL, y = "metabolites Shannon index") +
  theme_bw()
shannon_lineplot
# Richness violin plot
richness_violinplot <- ggviolin(richness_table, x = "Group", y = "richness",
                                color = "Group", fill = "Group", palette = c("grey", "#c97937", "#3b6291"), add = "mean_sd",
                                add.params = list(color = "white", size = 0.5)) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif", paired = TRUE, p.adjust.method1 = "BH") +
  labs(y = "metabolites richness", x = "") +
  theme_bw()
# Richness line plot
richness_lineplot <- ggline(richness_table, x = "Group", y = "richness", group = "Volunteer_ID",
                            size = 0.3, color = "Volunteer_ID", palette = "lancet") +
  labs(x = NULL, y = "metabolites richness") +
  theme_bw()
richness_lineplot
# Evenness violin plot
evenness_violinplot <- ggviolin(evenness_table, x = "Group", y = "evenness",
                                color = "Group", fill = "Group", palette = c("grey", "#c97937", "#3b6291"), add = "mean_sd",
                                add.params = list(color = "white", size = 0.5)) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif", paired = TRUE, p.adjust.method1 = "BH") +
  labs(y = "metabolites evenness", x = "") +
  theme_bw()
# Evenness line plot
evenness_lineplot <- ggline(evenness_table, x = "Group", y = "evenness", group = "Volunteer_ID",
                            size = 0.3, color = "Volunteer_ID", palette = "lancet") +
  labs(x = NULL, y = "metabolites evenness") +
  theme_bw()
evenness_lineplot
# Combined plots
alpha_diversity_violinplot <- ggarrange(
  ggarrange(shannon_violinplot, richness_violinplot, evenness_violinplot, 
            nrow = 3, ncol = 1, common.legend = T, legend = "right"),
  ggarrange(shannon_lineplot, richness_lineplot, evenness_lineplot, 
            nrow = 3, ncol = 1, common.legend = T, legend = "right"),
  nrow = 1, ncol = 2)
alpha_diversity_violinplot  # Save as: Metabolite α diversity

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
species_pcoa$Group <- factor(species_pcoa$Group, levels = c("Baseline", "Adverse_reaction_occurred", "Tolerance_established"), ordered = T)
species_pcoa_plot = ggscatter(species_pcoa, x = "PC1", y = "PC2", 
                                 color = "Group", palette = c("grey","#c97937","#3b6291"), ellipse = TRUE, ellipse.type = "norm",
                                 ellipse.border.remove = F, ellipse.alpha = 0.05, ellipse.level = 0.7, 
                                 size = 1,  legend = "right", mean.point =T, mean.point.size = 4,
                                 xlab = paste0("PCoA1 ( ",pc1,"%"," )",sep="" ),
                                 ylab = paste0("PCoA2 ( ",pc2,"%"," )",sep="" ))+theme_bw()
species_pcoa_plot   # Save as: Metabolite PCoA plot
######## adonis dissimilarity analysis
sink('results/metabolite_analysis/metabolites_dissimilarity_test_results.txt', split = T)
set.seed(1234)
pairwise.adonis2(data1~Group, data, nperm = 9999)
sink()


# --------------------------------------New analysis, selecting TOP500 metabolites from intolerant group and significant baseline differences between groups---------------------------------

# Baseline significance testing
    skin_trait <- read.csv("results/metabolite_analysis/metabolities_retinol_yes_no.csv", header=T, sep=",", check.names = F) #txt, tsv  #rows as samples, columns as factors
    dim(skin_trait)
    head(skin_trait)
    skin_trait$sample <- factor(skin_trait$sample)
    skin_trait$Group <- factor(skin_trait$Group, levels = c("Intolerant_Baseline", "Tolerant_Baseline"), ordered = T)
    colnames(skin_trait[,-c(1:2)])[1]
    #Group-wise statistical calculations
    skin_trait_melt_table <- reshape2::melt(skin_trait, id = 1:2, variable.name = "skin_trait_index", value.name = "Value")
    dim(skin_trait_melt_table)
    head(skin_trait_melt_table,20)
    #calculate mean_sd
    skin_trait_mean_sd_results <- skin_trait_melt_table %>% group_by(skin_trait_index, Group) %>% summarise(mean = mean(Value), sd = sd(Value))
    write.csv(skin_trait_mean_sd_results, "results/metabolite_analysis/compare_metabolites_abundance_mean_sd_Group.csv")
    #calculate significance
    skin_trait_stat_test <- compare_means(Value ~ Group, group.by = "skin_trait_index", 
                                          #paired = T,
                                          data = skin_trait_melt_table)
    write.csv(skin_trait_stat_test, "results/metabolite_analysis/compare_metabolites_stat_test.csv")

# For metabolites.csv, filter to retain metabolites with significant baseline differences and beneficial/harmful metabolites (already implicitly TOP500)

    # Beneficial metabolites: Sphinganine, SM(d18:0/20:2(11Z,14Z)), Malic acid, Aldosterone 18-glucuronide
    # Harmful metabolites: PC(22:4(7Z,10Z,13Z,16Z)/22:5(4Z,7Z,10Z,13Z,16Z)), Dihydroxyacetone Phosphate Acyl Ester, LysoPC(22:6(4Z,7Z,10Z,13Z,16Z,19Z)/0:0), 9,10-Epoxyoctadecenoic acid   

    stat_test <- read.csv("results/metabolite_analysis/compare_metabolites_stat_test.csv", header = TRUE, stringsAsFactors = FALSE)
    sig_threshold <- 0.05

    # Filter metabolites with significant differences
    sig_metabolites <- stat_test$skin_trait_index[stat_test$p < sig_threshold]
    cat("Number of significantly different metabolites:", length(sig_metabolites), "\n")

    # Define beneficial and harmful metabolites
    beneficial_metabolites <- c("Sphinganine", "SM(d18:0/20:2(11Z,14Z))", "Malic Acid", "Aldosterone 18-Glucuronide")
    harmful_metabolites <- c("PC(22:4(7Z,10Z,13Z,16Z)/22:5(4Z,7Z,10Z,13Z,16Z))", 
                            "Dihydroxyacetone Phosphate Acyl Ester", 
                            "PC(22:6(4Z,7Z,10Z,13Z,16Z,19Z)/0:0)", # LysoPC(22:6(4Z,7Z,10Z,13Z,16Z,19Z)/0:0) not found, used PC(22:6(4Z,7Z,10Z,13Z,16Z,19Z)/0:0) (also TOP500)
                            "(9s,10s)-9,10-Dihydroxyoctadecanoic Acid")  # 9,10-Epoxyoctadecenoic acid not found, used (9S,10S)-9,10-Dihydroxyoctadecanoic Acid (also TOP500, hydrolysis product of the former)

    # All beneficial and harmful metabolites
    interest_metabolites <- c(beneficial_metabolites, harmful_metabolites)
    cat("Number of beneficial and harmful metabolites of interest:", length(interest_metabolites), "\n")

    # Baseline significantly different and beneficial or harmful metabolites
    common_metabolites <- intersect(sig_metabolites, interest_metabolites)
    cat("Number of intersecting metabolites:", length(common_metabolites), "\n")
    cat("List of intersecting metabolites:\n")
    print(common_metabolites)

    selected_metabolites <- common_metabolites

    # Print
    if(length(selected_metabolites) > 0) {
      cat("List of metabolites meeting criteria (beneficial/harmful):\n")
      for(i in 1:length(selected_metabolites)) {
        if(selected_metabolites[i] %in% beneficial_metabolites) {
          cat(i, ". ", selected_metabolites[i], " (beneficial)\n", sep="")
        } else {
          cat(i, ". ", selected_metabolites[i], " (harmful)\n", sep="")
        }
      }
    } else {
      cat("No metabolites satisfy belonging to beneficial/harmful metabolite list\n")
    }

    data_met <- read.csv("metabolites.csv", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

    # Get the first 3 columns as sample information columns
    sample_info_cols <- 1:3

    # First create a dataframe containing only sample information columns
    filtered_data <- data_met[, sample_info_cols]

    # Add metabolite columns in the order of original definition
    found_metabolites <- c()
    for(metabolite in interest_metabolites) {
      if(metabolite %in% colnames(data_met)) {
        # Add found metabolites to new dataframe
        filtered_data[[metabolite]] <- data_met[[metabolite]]
        found_metabolites <- c(found_metabolites, metabolite)
      }
    }

    cat("\nNumber of beneficial or harmful metabolites found in dataset:", length(found_metabolites), "\n")

    if(length(found_metabolites) > 0) {
      # List found beneficial or harmful metabolites (in original order)
      cat("Beneficial or harmful metabolites found in dataset:\n")
      for(i in 1:length(found_metabolites)) {
        met_name <- found_metabolites[i]
        if(met_name %in% beneficial_metabolites) {
          cat(i, ". ", met_name, " (beneficial)\n", sep="")
        } else {
          cat(i, ". ", met_name, " (harmful)\n", sep="")
        }
      }
      
      # Save results to new file
      write.csv(filtered_data, "results/metabolite_analysis/filtered_beneficial_harmful_metabolites.csv", row.names = FALSE)
      cat("\nSuccessfully filtered and saved beneficial or harmful metabolite data (in original definition order)\n")
    }

# Plot metabolites meeting the above conditions
    metabolites_data <- read.csv("results/metabolite_analysis/filtered_beneficial_harmful_metabolites.csv", 
                            header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    # Set factor variables
    metabolites_data$sample <- factor(metabolites_data$sample)
    metabolites_data$Group <- factor(metabolites_data$Group, 
                                  levels = c("Baseline", "Adverse_reaction_occurred", "Tolerance_established"), 
                                  ordered = TRUE)
    metabolites_data$Volunteer_ID <- factor(metabolites_data$Volunteer_ID)

    # Convert data to long format
    library(reshape2)
    metabolites_melt <- reshape2::melt(metabolites_data, 
                                      id = 1:3, 
                                      variable.name = "skin_trait_index", 
                                      value.name = "Value")

    # Set comparison groups
    comparisons <- list(c("Baseline", "Adverse_reaction_occurred"), 
                        c("Adverse_reaction_occurred", "Tolerance_established"),
                        c("Baseline", "Tolerance_established"))

    skin_trait_violinplot <- ggboxplot(
      metabolites_melt, 
      x = "Group", 
      y = "Value",
      color = "Group", 
      fill = "Group", 
      palette = c("grey", "#c97937", "#3b6291"), 
      add = "mean_sd",
      add.params = list(color = "black", size = 0.2)
    ) +
      theme_bw() +
      labs(x = "", y = "") +
      facet_wrap(
        ~skin_trait_index, 
        ncol = 4, 
        nrow = 2, 
        scales = "free_y"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.8),
        strip.text = element_text(size = 8)
      ) +
      stat_compare_means(
        method = "wilcox.test", 
        paired = FALSE, 
        label = "p.signif", 
        hide.ns = TRUE, 
        label.y.npc = "top", 
        vjust = 0.5, 
        comparisons = comparisons, 
        p.adjust.method = "BH"
      )

    print(skin_trait_violinplot)  # Baseline significantly different and beneficial/harmful metabolites


# Significance testing of these 8 metabolites in tolerance group (4 appeared in tolerance group)
# Read new metabolite data
metabolites_data <- read.csv("results/metabolite_analysis/metabolites_8meta.csv", 
                            header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
# Check original timepoints
print("Original unique values in Time column:")
print(unique(metabolites_data$Time))
# Keep only needed timepoints
metabolites_data <- metabolites_data[metabolites_data$Time %in% c("Day0", "Day7", "Day21"), ]
# Check timepoints after filtering
print("Unique values in Time column after filtering:")
print(unique(metabolites_data$Time))
# Set factor variables
metabolites_data$sample <- factor(metabolites_data$sample)
metabolites_data$Time <- factor(metabolites_data$Time, 
                               levels = c("Day0", "Day7", "Day21"), 
                               ordered = TRUE)
metabolites_data$Volunteer_ID <- factor(metabolites_data$Volunteer_ID)
# Convert data to long format
library(reshape2)
metabolites_melt <- reshape2::melt(metabolites_data, 
                                  id = 1:3,  # sample, Time, Volunteer_ID
                                  variable.name = "skin_trait_index", 
                                  value.name = "Value")
# Set comparison groups (changed to timepoint comparisons)
comparisons <- list(c("Day0", "Day7"), 
                    c("Day7", "Day21"),
                    c("Day0", "Day21"))
# Create boxplot
skin_trait_violinplot <- ggboxplot(
  metabolites_melt, 
  x = "Time",  # Changed to Time
  y = "Value",
  color = "Time",  # Changed to Time
  fill = "Time",   # Changed to Time
  palette = c("#2E8B57", "#FF6347", "#4682B4"),  # 3 colors for 3 timepoints
  add = "mean_sd",
  add.params = list(color = "black", size = 0.2)
) +
  theme_bw() +
  labs(x = "", y = "") +
  facet_wrap(
    ~skin_trait_index, 
    ncol = 2,  # Adjust according to your 4 metabolites
    nrow = 2, 
    scales = "free_y"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.8),
    strip.text = element_text(size = 8)
  ) +
  stat_compare_means(
    method = "wilcox.test", 
    paired = FALSE,  # Can change to TRUE if paired data
    label = "p.signif", 
    hide.ns = FALSE, 
    vjust = 0.3,  # Adjust label position
    step.increase = 0.1,  # Adjust spacing between multiple comparisons
    comparisons = comparisons, 
    p.adjust.method = "BH",
    size = 3
  )
print(skin_trait_violinplot)  # Significance testing of 8 tolerance group metabolites at three timepoints