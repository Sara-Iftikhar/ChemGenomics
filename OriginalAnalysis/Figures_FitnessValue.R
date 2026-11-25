# -----------------------------------------------------------------------------
# Script: Fitness value figures and downstream exploratory analyses
# Purpose: Generate figures for colony-size/opacity relationships, replicate
#          consistency, condition categories, GWAS summaries, plasmid effects,
#          machine-learning performance, and ST/virulence/resistance analyses.
# Data root: all project paths have been updated to use /Data/.
# Note: code logic is unchanged; only comments and requested path updates were added.
# -----------------------------------------------------------------------------
#Results scatter plot Figure 1
# Load core tidyverse functions for data import, wrangling, and plotting.
library(tidyverse)
# Directory containing IRIS colony-measurement output files.
dir_path <- "/Data/assay_results_full"
files <- list.files(dir_path,pattern = "-1_A\\.JPG\\.iris$")
files<-gsub("-1_A.JPG.iris","",files)
conditions<- unique(as.character( sapply(files, function(x) paste0( strsplit(x,"-")[[1]][1], "-",  strsplit(x,"-")[[1]][2]) ) ))  

# Containers for replicate-level size-opacity correlations and pooled observations.
cor_size_opacity_list <- c()
cor_size_total_df <- c()
# Loop through each condition and read the four plate/replicate IRIS files.
for(i in 1:length(conditions)){
  print(i/length(conditions))
  condition_input<-conditions[i]
  # Construct input paths for the four replicate/plate files for the current condition.
  m1_tag <- paste0("/Data/assay_results_full/", condition_input, "-1-1_A.JPG.iris")
  m2_tag <- paste0("/Data/assay_results_full/", condition_input, "-2-1_A.JPG.iris")
  m3_tag <- paste0("/Data/assay_results_full/", condition_input, "-3-1_A.JPG.iris")
  m4_tag <- paste0("/Data/assay_results_full/", condition_input, "-4-1_A.JPG.iris")
  
  # Read IRIS outputs while ignoring metadata/header comment lines.
  m1 <- read_tsv(m1_tag, comment = "#", show_col_types = FALSE)
  m2 <- read_tsv(m2_tag, comment = "#", show_col_types = FALSE)
  m3 <- read_tsv(m3_tag, comment = "#", show_col_types = FALSE)
  m4 <- read_tsv(m4_tag, comment = "#", show_col_types = FALSE)
  
  # Calculate colony-size versus opacity correlation for each replicate.
  cor_size_opacity_m1 <- cor(m1$`colony size`, m1$opacity)
  cor_size_opacity_m2 <- cor(m2$`colony size`, m2$opacity)
  cor_size_opacity_m3 <- cor(m3$`colony size`, m3$opacity)
  cor_size_opacity_m4 <- cor(m4$`colony size`, m4$opacity)
  
  # Pool colony-size and opacity values across all four replicates.
  tmp_size <- c(m1$`colony size`, m2$`colony size`,m3$`colony size`, m4$`colony size`)
  tmp_opacity <- c(m1$opacity, m2$opacity,m3$opacity, m4$opacity)
  
  cor_size_total_df <- rbind(cor_size_total_df, cbind(tmp_size,tmp_opacity))
  
  cor_size_opacity_list <- c(cor_size_opacity_list,cor_size_opacity_m1,cor_size_opacity_m2,cor_size_opacity_m3,cor_size_opacity_m4)
}
# Summarise replicate-level size-opacity correlations across all conditions.
mean(cor_size_opacity_list)

# Convert pooled observations to a data frame for plotting.
cor_size_total_df <- data.frame(cor_size_total_df)
colnames(cor_size_total_df)<- c("Size","Opacity")
head(cor_size_total_df)

# Randomly down-sample observations to make density plotting lighter.
cor_size_total_df_short <- cor_size_total_df[
  sample(nrow(cor_size_total_df), size = round(0.1 * nrow(cor_size_total_df))),
]

# Load density-aware scatter plotting functions.
library(ggpointdensity)

cor_size_total_df_short$Size <- cor_size_total_df_short$Size/max(cor_size_total_df_short$Size)
cor_size_total_df_short$Opacity <- cor_size_total_df_short$Opacity/max(cor_size_total_df_short$Opacity)

# Plot relative colony size versus opacity, colored by local point density.
ggplot(cor_size_total_df_short , aes(Size, Opacity)) +
  geom_pointdensity(
    aes(alpha = sqrt(after_stat(density))),
    size = 0.05, method = "neighbors"
  ) +
  scale_color_viridis_c(option = "inferno") +
  scale_alpha(range = c(0.05, 1), guide = "none") +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed",
              color = "red",
              linewidth = 0.5) +
  ylab("Rel. Opacity")+
  xlab("Rel. Size")+
 # coord_equal(xlim = lims, ylim = lims, expand = FALSE) +
  theme_bw() +
  ylim(range(0,1))+
  xlim(range(0,1))+
  theme(
    axis.text.x = element_text(size = 15, angle = 0, hjust = 1),
    axis.text.y = element_text(size = 15, hjust = 1),
    axis.title.x = element_text(color = "black", size = 18, face = "bold"),
    axis.title.y = element_text(color = "black", size = 18, face = "bold"),
    strip.text.x = element_text(size = 17, color = "black", face = "bold")
  )


# -----------------------------------------------------------------------------
# Figure 2: replicate-level colony-size consistency across plates.
# This section links IRIS measurements to plate-layout metadata, reshapes
# repeated strain measurements into wide format, and compares replicates.
# -----------------------------------------------------------------------------
#Figure 2
# Accumulate strain-level replicate measurements across all conditions.
m_tot_df<-c()
for(i in 1:length(conditions)){
  #for(i in 1:10){
  
  print(i/length(conditions))
  condition_input<-conditions[i]
  # Read plate-layout metadata that maps row/column positions to strain IDs.
  p1 <- read_tsv("/Data/conditions/plate_tr_1.tsv", show_col_types = FALSE)
  p2 <- read_tsv("/Data/conditions/plate_tr_2.tsv", show_col_types = FALSE)
  p3 <- read_tsv("/Data/conditions/plate_tr_3.tsv", show_col_types = FALSE)
  p4 <- read_tsv("/Data/conditions/plate_tr_4.tsv", show_col_types = FALSE)
  
  m1_tag <- paste0("/Data/assay_results_full/", condition_input, "-1-1_A.JPG.iris")
  m2_tag <- paste0("/Data/assay_results_full/", condition_input, "-2-1_A.JPG.iris")
  m3_tag <- paste0("/Data/assay_results_full/", condition_input, "-3-1_A.JPG.iris")
  m4_tag <- paste0("/Data/assay_results_full/", condition_input, "-4-1_A.JPG.iris")
  
  m1 <- read_tsv(m1_tag, comment = "#", show_col_types = FALSE)
  m2 <- read_tsv(m2_tag, comment = "#", show_col_types = FALSE)
  m3 <- read_tsv(m3_tag, comment = "#", show_col_types = FALSE)
  m4 <- read_tsv(m4_tag, comment = "#", show_col_types = FALSE)
  
  # Standardize column names so they can be referenced without spaces.
  colnames(m1) <- gsub(" ", "_", colnames(m1))
  colnames(m2) <- gsub(" ", "_", colnames(m2))
  colnames(m3) <- gsub(" ", "_", colnames(m3))
  colnames(m4) <- gsub(" ", "_", colnames(m4))
  
  # Create row-column keys for joining IRIS measurements to plate metadata.
  m1$r_c <- paste(m1$row, m1$column)
  m2$r_c <- paste(m2$row, m2$column)
  m3$r_c <- paste(m3$row, m3$column)
  m4$r_c <- paste(m4$row, m4$column)
  
  # Attach strain labels to each colony measurement using row-column keys.
  m1 <- left_join(m1, p1, by = c("r_c" = "r_c"))
  m2 <- left_join(m2, p2, by = c("r_c" = "r_c"))
  m3 <- left_join(m3, p3, by = c("r_c" = "r_c"))
  m4 <- left_join(m4, p4, by = c("r_c" = "r_c"))
  
  # Keep only positional information, colony size, and strain identity.
  m1 <- dplyr::select(m1, row.x, column.x, colony_size, strain)
  m2 <- dplyr::select(m2, row.x, column.x, colony_size, strain)
  m3 <- dplyr::select(m3, row.x, column.x, colony_size, strain)
  m4 <- dplyr::select(m4, row.x, column.x, colony_size, strain)
  
  m1_zero <- m1[m1$colony_size == 0, ]
  m2_zero <- m2[m2$colony_size == 0, ]
  m3_zero <- m3[m3$colony_size == 0, ]
  m4_zero <- m4[m4$colony_size == 0, ]
  
  # Retain positive colony-size measurements for replicate comparison.
  m1 <- m1[m1$colony_size > 0, ]
  m2 <- m2[m2$colony_size > 0, ]
  m3 <- m3[m3$colony_size > 0, ]
  m4 <- m4[m4$colony_size > 0, ]
  
  m1 <- m1[,-c(1,2)]
  # Reshape replicate colony-size measurements into one row per strain.
  m1_wide <- m1 %>%
    group_by(strain) %>%
    mutate(meas_id = row_number()) %>%   # create measurement index
    ungroup() %>%
    pivot_wider(
      names_from  = meas_id,
      values_from = colony_size,
      names_prefix = "meas"
    )
  m1_wide <- m1_wide[,1:5]
  
  m2 <- m2[,-c(1,2)]
  m2_wide <- m2 %>%
    group_by(strain) %>%
    mutate(meas_id = row_number()) %>%   # create measurement index
    ungroup() %>%
    pivot_wider(
      names_from  = meas_id,
      values_from = colony_size,
      names_prefix = "meas"
    )
  m2_wide <- m2_wide[,1:5]
  
  
  m3 <- m3[,-c(1,2)]
  m3_wide <- m3 %>%
    group_by(strain) %>%
    mutate(meas_id = row_number()) %>%   # create measurement index
    ungroup() %>%
    pivot_wider(
      names_from  = meas_id,
      values_from = colony_size,
      names_prefix = "meas"
    )
  m3_wide <- m3_wide[,1:5]
  
  
  m4 <- m4[,-c(1,2)]
  m4_wide <- m4 %>%
    group_by(strain) %>%
    mutate(meas_id = row_number()) %>%   # create measurement index
    ungroup() %>%
    pivot_wider(
      names_from  = meas_id,
      values_from = colony_size,
      names_prefix = "meas"
    )
  m4_wide <- m4_wide[,1:5]
  
  # Combine all plates for this condition before appending to the global table.
  m_tot<- rbind(m1_wide,m2_wide,m3_wide,m4_wide)
  m_tot_df<-rbind(m_tot_df,m_tot)
}

# Keep only strains with complete replicate measurements.
m_tot_df_short <- m_tot_df[complete.cases(m_tot_df),]


library(ggpointdensity)
# Quick replicate comparison plot for the first pair of measurements.
ggplot(m_tot_df_short, aes(meas1, meas2)) +
  geom_pointdensity(
    aes(alpha = sqrt(after_stat(density))),
    size = 0.4   # ← smaller points
  ) +
  scale_color_viridis_c(option = "inferno") +
  scale_alpha(range = c(0.05, 1), guide = "none") +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed",
              color = "grey",
              linewidth = 0.7) +
  ylab("Rel. Opacity")+
  xlab("Rel. Size")+
  coord_equal() +
  theme_classic()


# Compute common axis limits so replicate plots are directly comparable.
# Compute common limits
lims <- range(c(m_tot_df_short$meas1,
                m_tot_df_short$meas2),
              na.rm = TRUE)

lims <- range(c(0,
                10000),
              na.rm = TRUE)

# Publication-style density plot comparing replicate measurements 3 and 4.
ggplot(m_tot_df_short, aes(meas3, meas4)) +
  geom_pointdensity(
    aes(alpha = sqrt(after_stat(density))),
    size = 0.01, method = "neighbors"
  ) +
  scale_color_viridis_c(option = "inferno") +
  scale_alpha(range = c(0.025, 1), guide = "none") +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed",
              color = "red",
              linewidth = 0.5) +
  coord_equal(xlim = lims, ylim = lims, expand = FALSE) +
  ylab("Replicate 1")+
  xlab("Replicate 2")+
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 15, angle = 0, hjust = 1),
    axis.text.y = element_text(size = 15, hjust = 1),
    axis.title.x = element_text(color = "black", size = 18, face = "bold"),
    axis.title.y = element_text(color = "black", size = 18, face = "bold"),
    strip.text.x = element_text(size = 17, color = "black", face = "bold")
  )



# Use numeric replicate columns to estimate pairwise replicate correlations.
# select only numeric measurement columns
meas_mat <- m_tot_df[, c("meas1", "meas2", "meas3", "meas4")]

# Spearman correlation assesses monotonic agreement between replicate measurements.
cor_matrix <- cor(meas_mat, method = "spearman", use = "pairwise.complete.obs")

cor_matrix[lower.tri(cor_matrix, diag = F)]


# -----------------------------------------------------------------------------
# Transition into condition-category summaries and ICC-related exploratory notes.
# -----------------------------------------------------------------------------
#Figure ICC distribution per codntion and aggregate 

#589672 SAR
#191887 SAR
#563000 SAR
132,120.15 USD

#30865


# -----------------------------------------------------------------------------
# Condition-category composition plots.
# This section aligns condition metadata to fitness-matrix columns and visualizes
# broad categories and chemical-class subgroups.
# -----------------------------------------------------------------------------
#Pie chart 
fintess_mat <- read_csv("/Data/absolute_fitness_ML.csv")
colnames(fintess_mat)


#CeferdericolCeftazadime-0125,1ugml
#CeferdericolTigecycline-0625,2ugml
#ColistinMeropenem-02,006ugml
#ColistinMeropenem-04,006ugml

fintess_mat <- read_csv("/Data/absolute_fitness_ML.csv")
colnames(fintess_mat)

condition_group<- read_csv("/Data/condition_group_11March.csv",show_col_types = FALSE)
condition_group$File_name

# Clean file-name labels so they match sanitized fitness-matrix column names.
condition_group$File_name <- gsub("[^A-Za-z0-9_]", "", condition_group$File_name )
condition_group <- condition_group[match(colnames(fintess_mat)[-1],condition_group$File_name),]

#condition_group <- condition_group[condition_group$Category!="Base",]
condition_group <- condition_group[!is.na(condition_group$Chemical),]

colnames(condition_group)
condition_group$Chemical_class
head(condition_group)


library(ggplot2)

# Count frequency of each Category
# Count how many conditions belong to each broad category.
category_counts <- condition_group %>%
  count(Category) %>%
  mutate(percentage = n / sum(n) * 100)

# Consistent condition-category color palette reused across several figures.
category_colors <- c(
  "Chemical stress" = "#377EB8",
  "Carbon source utilisation" = "#4DAF4A",
  "Antibiotics" = "#E41A1C",
  "Environmental stress" = "#984EA3",
  "Antiseptic" = "#FF7F00",
  "Stress combination" = "#A65628",
  "Metal stress" = "#B8860B",
  "Base" = "#999999",
  "Envelope stress" = "#17BECF"
)

# Pie chart showing counts by broad condition category.
ggplot(category_counts, aes(x = "", y = n, fill = Category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  
  scale_fill_manual(values = category_colors) +  # 🔴 this is the key line
  geom_text(
    aes(label = n),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "black"
  )+
  labs(
       fill = "Category") +
  
  theme_void() +
  
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )



library(dplyr)
library(ggplot2)
library(scales)

# Count categories
category_counts <- condition_group %>%
  count(Category) %>%
  mutate(
    percentage = n / sum(n),
    label = paste0(Category, "\n", percent(percentage, accuracy = 0.1))
  )

# Donut chart
# Donut chart with percentages
# Donut chart version of the category count visualization.
ggplot(category_counts, aes(x = 2, y = n, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +   # controls hole size
  geom_text(
    aes(label = n),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "black"
  ) +
  scale_fill_manual(values = category_colors) +  # 🔴 this is the key line
  theme_void() +
  labs(
    title = "",
    fill = "Category"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )




ggplot(category_counts, aes(x = 2, y = n, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +   # controls hole size
  geom_text(
    aes(label = scales::percent(n / sum(n), accuracy = 0.1)),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "black"
  ) +
  scale_fill_manual(values = c(
    "#E41A1C",  # red
    "#377EB8",  # blue
    "#4DAF4A",  # green
    "#984EA3",  # purple
    "#FF7F00",  # orange
    "#A65628",  # brown
    "#F781BF",  # pink
    "#999999",   # gray
    "#00FFFF"   # gray
  )) +
  theme_void() +
  labs(
    title = "Category Distribution",
    fill = "Category"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )




# -----------------------------------------------------------------------------
# Chemical stress subgrouping: chemical classes are collapsed into broader
# interpretable mechanisms before plotting.
# -----------------------------------------------------------------------------
#stress

#stresses
#fintess_mat <- read_csv("/Data/absolute_fitness_ML.csv")
#colnames(fintess_mat)

#condition_group<- read_csv("/Data/condition_group.csv",show_col_types = FALSE)
#condition_group$File_name

#condition_group$File_name <- gsub("[^A-Za-z0-9_]", "", condition_group$File_name )
#condition_group <- condition_group[match(colnames(fintess_mat)[-1],condition_group$File_name),]

# Subset chemical-stress conditions for subgrouping.
condition_group_short<- condition_group[condition_group$Category=="Chemical stress",]
table(condition_group_short$Chemical_class)
condition_group_short$Chemical

condition_group_short <- condition_group_short %>%
  mutate(Chemical_class = case_when(
    
    Chemical_class %in% c(
      "Alcohol",
      "Aliphatic amine",
      "Aromatic heterocycle (tryptophan catabolite)",
      "Aromatic phenol",
      "Flavonoid (polyphenol)",
      "Organic acid",
      "Polar aprotic solvent",
      "Condition (acidic environment)"
    ) ~ "Organic small molecule",
    
    Chemical_class %in% c(
      "Anionic surfactant",
      "Bile acid salt",
      "Ionophore (uncoupler)"
    ) ~ "Surfactant / membrane-active",
    
    Chemical_class %in% c(
      "Heavy-metal salt"
    ) ~ "Metal / inorganic stress",
    
    Chemical_class %in% c(
      "Biguanide",
      "Polypeptide (antimicrobial)",
      "Thiourea derivative",
      "Phenylalkylamine"
    ) ~ "Antimicrobial / bioactive",
    
    TRUE ~ "Other"
  ))



# Count categories
category_counts <- condition_group_short %>%
  count(Chemical_class) %>%
  mutate(
    percentage = n / sum(n),
    label = paste0(Chemical_class, "\n", percent(percentage, accuracy = 0.1))
  )

# Donut chart
ggplot(category_counts, aes(x = 2, y = n, fill = Chemical_class)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +   # controls hole size (increase 0.5 to enlarge hole)
  theme_void() +
  labs(title = "",
       fill = "Category") +
  geom_text(
    aes(label =n),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "black"
  ) +
  scale_fill_manual(values = c(
    "#E41A1C",  # red
    "#377EB8",  # blue
    "#4DAF4A",  # green
    "#999999",   # gray
    "#984EA3",  # purple
    
    "#00FFFF"   # gray
  )) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )

#stresses
#condition_group<- read_csv("/Data/condition_group.csv")

#condition_group$File_name <- gsub("[^A-Za-z0-9_]", "", condition_group$File_name )
#condition_group <- condition_group[match(colnames(fintess_mat)[-1],condition_group$File_name),]

# Subset carbon-source-utilisation conditions for subgrouping.
condition_group_short<- condition_group[condition_group$Category=="Carbon source utilisation",]
table(condition_group_short$Chemical_class)
condition_group_short$Chemical
condition_group_short <- condition_group_short %>%
  mutate(Chemical_class = case_when(
    
    Chemical_class %in% c(
      "Carbohydrate",
      "Disaccharide",
      "Disaccharide (carbohydrate)",
      "Monosaccharide",
      "Sugar alcohol"
    ) ~ "Carbon source",
    
    Chemical_class %in% c(
      "Amino acid",
      "Nucleoside (pyrimidine)"
    ) ~ "Nitrogen / nucleic substrate",
    
    Chemical_class %in% c(
      "Biological fluid",
      "Biological fluid (growth medium)",
      "Complex serum (growth medium)"
    ) ~ "Complex nutrient environment",
    
    Chemical_class %in% c(
      "Organic chaotrope (urea)",
      "Osmotic condition"
    ) ~ "Stress / physicochemical condition",
    
    TRUE ~ "Other"
  ))


# Count categories
category_counts <- condition_group_short %>%
  count(Chemical_class) %>%
  mutate(
    percentage = n / sum(n),
    label = paste0(Chemical_class, "\n", percent(percentage, accuracy = 0.1))
  )

# Donut chart
ggplot(category_counts, aes(x = 2, y = n, fill = Chemical_class)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +   # controls hole size (increase 0.5 to enlarge hole)
  theme_void() +
  labs(title = "",
       fill = "Category") +
  geom_text(
    aes(label =n),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "black"
  ) +
  scale_fill_manual(values = c(
    "#E41A1C",  # red
    "#377EB8",  # blue
    "#4DAF4A",  # green
    "#999999",   # gray
    "#984EA3",  # purple
    
    "#00FFFF"   # gray
  )) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )



# -----------------------------------------------------------------------------
# Environmental stress subgrouping: conditions are grouped into pH, ionic/salt,
# osmotic, and physical/environmental classes.
# -----------------------------------------------------------------------------
#Environmental stress
#condition_group<- read_csv("/Data/condition_group.csv")

#condition_group$File_name <- gsub("[^A-Za-z0-9_]", "", condition_group$File_name )
#condition_group <- condition_group[match(colnames(fintess_mat)[-1],condition_group$File_name),]

# Subset environmental-stress conditions for subgrouping.
condition_group_short<- condition_group[condition_group$Category=="Environmental stress",]
table(condition_group_short$Chemical_class)
condition_group_short$Chemical


condition_group_short <- condition_group_short %>%
  mutate(Chemical_class = case_when(
    
    Chemical_class %in% c(
      "Buffer (morpholine propanesulfonic acid)",
      "Condition (acidic environment)",
      "Condition (alkaline environment)",
      "Condition (near-neutral pH)",
      "Condition (strongly alkaline)"
    ) ~ "pH / buffering condition",
    
    Chemical_class %in% c(
      "Inorganic salt (carbonate)",
      "Inorganic salt (nitrate)",
      "Ionic compound (salt)"
    ) ~ "Ionic / salt stress",
    
    Chemical_class %in% c(
      "Osmotic condition"
    ) ~ "Osmotic stress",
    
    Chemical_class %in% c(
      "Condition (no oxygen)",
      "Physical condition",
      "Radiation"
    ) ~ "Physical / environmental stress",
    
    TRUE ~ "Other"
  ))


# Count categories
category_counts <- condition_group_short %>%
  count(Chemical_class) %>%
  mutate(
    percentage = n / sum(n),
    label = paste0(Chemical_class, "\n", percent(percentage, accuracy = 0.1))
  )

# Donut chart
ggplot(category_counts, aes(x = 2, y = n, fill = Chemical_class)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +   # controls hole size (increase 0.5 to enlarge hole)
  theme_void() +
  labs(title = "",
       fill = "Category") +
  geom_text(
    aes(label =n),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "black"
  ) +
  scale_fill_manual(values = c(
    "#E41A1C",  # red
    "#377EB8",  # blue
    "#999999",   # gray
    "#4DAF4A",  # green
    "#984EA3",  # purple
    "#00FFFF"   # gray
  )) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )


# -----------------------------------------------------------------------------
# Stress-combination subgrouping: classify combinations by antibiotic,
# non-antibiotic, and environmental-condition components.
# -----------------------------------------------------------------------------
#stress combination 
#condition_group<- read_csv("/Data/condition_group.csv")
#sort(table(condition_group$Category))

#condition_group$File_name <- gsub("[^A-Za-z0-9_]", "", condition_group$File_name )
#condition_group <- condition_group[match(colnames(fintess_mat)[-1],condition_group$File_name),]

# Subset stress-combination conditions for subgrouping.
condition_group_short<- condition_group[condition_group$Category=="Stress combination",]
table(condition_group_short$Chemical)

condition_group_short <- condition_group_short %>%
  mutate(Chemical_class = case_when(
    
    Chemical %in% c(
      "CefepimeTazobactam",
      "CefiderocolCeftazadime",
      "CefiderocolColistin",
      "CefiderocolImipenem",
      "CefiderocolMeropenem",
      "CefiderocolTigecycline",
      "CiprofloxacinRifampicin",
      "ColistinCefepime",
      "ColistinCiprofloxacin",
      "ColistinGentamicin",
      "ColistinMeropenem",
      "PiperacillinTazobactam"
    ) ~ "Antibiotic–antibiotic combination",
    
    Chemical %in% c(
      "ColistinMetformin"
    ) ~ "Antibiotic–non-antibiotic combination",
    
    Chemical %in% c(
      "ColistinpH5",
      "ColistinpH9",
      "TigecyclinepH5"
    ) ~ "Antibiotic–environment condition",
    
    TRUE ~ "Other"
  ))

category_counts <- condition_group_short %>%
  count(Chemical_class) %>%
  mutate(
    percentage = n / sum(n),
    label = paste0(Chemical_class, "\n", percent(percentage, accuracy = 0.1))
  )


ggplot(category_counts, aes(x = 2, y = n, fill = Chemical_class)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +   # controls hole size (increase 0.5 to enlarge hole)
  theme_void() +
  labs(title = "",
       fill = "Category") +
  geom_text(
    aes(label =n),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "black"
  ) +
  scale_fill_manual(values = c(
    "#E41A1C",  # red
    "#377EB8",  # blue
    "#4DAF4A",  # green
    "#999999",   # gray

    "#984EA3",  # purple
    "#00FFFF"   # gray
  )) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )


# -----------------------------------------------------------------------------
# Antibiotic subgrouping: individual antibiotics are mapped to pharmacological
# classes for composition plots.
# -----------------------------------------------------------------------------
#antibiotcs
#condition_group<- read_csv("/Data/condition_group.csv")
#sort(table(condition_group$Category))

#condition_group$File_name <- gsub("[^A-Za-z0-9_]", "", condition_group$File_name )
#condition_group <- condition_group[match(colnames(fintess_mat)[-1],condition_group$File_name),]

# Subset antibiotic conditions for antibiotic-class grouping.
condition_group_short<- condition_group[condition_group$Category=="Antibiotics",]
table(condition_group_short$Chemical)

condition_group_short <- condition_group_short %>%
  mutate(Chemical_class = case_when(
    
    # Cephalosporins
    Chemical %in% c(
      "Cefaclor",
      "Cefepime",
      "Cefotaxime",
      "Ceftazidime",
      "Ceftriaxone",
      "Cephalexin"
    ) ~ "Cephalosporin",
    
    # Carbapenems
    Chemical %in% c(
      "Doripenem",
      "Ertapenem",
      "Imipenem",
      "Meropenem"
    ) ~ "Carbapenem",
    
    # Penicillins
    Chemical %in% c(
      "Ampicillin",
      "Oxacillin",
      "Piperacillin"
    ) ~ "Penicillin",
    
    # Monobactam-like siderophore cephalosporin
    Chemical %in% c(
      "Cefiderocol"
    ) ~ "Siderophore cephalosporin",
    
    # β-lactamase inhibitor
    Chemical %in% c(
      "Tazobactam"
    ) ~ "β-lactamase inhibitor",
    
    # Other antibiotic classes
    Chemical %in% c("Gentamicin", "Kanamycin") ~ "Aminoglycoside",
    Chemical %in% c("Ciprofloxacin") ~ "Fluoroquinolone",
    Chemical %in% c("Clindamycin") ~ "Lincosamide",
    Chemical %in% c("Chloramphenicol") ~ "Amphenicol",
    Chemical %in% c("Colistin (Polymyxin E)") ~ "Polymyxin",
    Chemical %in% c("Daptomycin") ~ "Lipopeptide",
    Chemical %in% c("Minocycline", "Tetracycline", "Tigecycline") ~ "Tetracycline",
    Chemical %in% c("Rifampicin") ~ "Rifamycin",
    Chemical %in% c("Trimethoprim") ~ "Folate pathway inhibitor",
    Chemical %in% c("Vancomycin") ~ "Glycopeptide",
    
    TRUE ~ "Other"
  ))



category_counts <- condition_group_short %>%
  count(Chemical_class) %>%
  mutate(
    percentage = n / sum(n),
    label = paste0(Chemical_class, "\n", percent(percentage, accuracy = 0.1))
  )

category_counts$Chemical_class <- factor(
  category_counts$Chemical_class,
  levels = c(
    
    # β-lactams grouped together
    "Penicillin",
    "Cephalosporin",
    "Carbapenem",
    "Siderophore cephalosporin",
    "β-lactamase inhibitor",
    
    # other classes
    "Aminoglycoside",
    "Fluoroquinolone",
    "Polymyxin",
    "Tetracycline",
    "Rifamycin",
    "Glycopeptide",
    "Lincosamide",
    "Lipopeptide",
    "Amphenicol",
    "Folate pathway inhibitor"
  )
)
ggplot(category_counts, aes(x = 2, y = n, fill = Chemical_class)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +   # controls hole size (increase 0.5 to enlarge hole)
  theme_void() +
  labs(title = "",
       fill = "Category") +
  geom_text(
    aes(label =n),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "black"
  ) +
  scale_fill_manual(values = c(
    
    # β-lactam subclasses (same blue family)
    "Cephalosporin" = "#4F81BD",
    "Carbapenem" = "#1F4E79",
    "Penicillin" = "#6FA8DC",
    "Siderophore cephalosporin"="#1F4EFF",
    "β-lactamase inhibitor" = "#9DC3E6",
    
    # Other antibiotic classes
    "Aminoglycoside" = "#4DAF4A",
    "Fluoroquinolone" = "#FF7F00",
    "Polymyxin" = "#E41A1C",
    "Tetracycline" = "#984EA3",
    "Rifamycin" = "#A65628",
    "Glycopeptide" = "#F781BF",
    "Folate pathway inhibitor" = "#999999",
    "Amphenicol" = "#66C2A5",
    "Lincosamide" = "#FFD92F",
    "Lipopeptide" = "#A6D854"
  )) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )


# -----------------------------------------------------------------------------
# Bimodality analysis.
# Hartigan dip statistics and density peak counts are calculated for each
# condition-specific median fitness distribution.
# -----------------------------------------------------------------------------
#Bimodality 
library(diptest)
# Directory containing IRIS colony-measurement output files.
dir_path <- "/Data/assay_results_full"
files <- list.files(dir_path,pattern = "-1_A\\.JPG\\.iris$")
files<-gsub("-1_A.JPG.iris","",files)

conditions<- unique(as.character( sapply(files, function(x) paste0( strsplit(x,"-")[[1]][1], "-",  strsplit(x,"-")[[1]][2]) ) ))  
# Containers for condition-level fitness distributions and modality summaries.
cond_df<-c()
d_dist<-c()
p_dist<-c()
sd_dist<-c()

for(i in 1:length(conditions)){
  print(i)
  cond<-read_csv(paste0(dir_path, "/", "condition_ICC_New_",conditions[i], ".csv"),show_col_types = FALSE)
  cond$condition<-conditions[i]
  cond_df<-rbind(cond_df, cond)
  
  # Store Hartigan dip statistic and p-value for multimodality assessment.
  d_dist<-c(d_dist,dip.test(cond$median_fc_overall)$statistic )
  p_dist<-c(p_dist,dip.test(cond$median_fc_overall)$p.value )
  
  h <- hist(cond$median_fc_overall, plot = FALSE)
  p <- h$counts / sum(h$counts)
  entropy <- -sum(p * log(p + 1e-12))  # Shannon entropy
  entropy
  #sd_dist<-c(sd_dist, entropy )
  
  f <- cond$median_fc_overall
  acf1 <- acf(f, plot = FALSE)$acf[2]
  acf1  # closer to 0 implies more ruggedness
  
  #sd_dist<-c(sd_dist, acf1)
  
  
  # Estimate density and count local peaks as an additional modality descriptor.
  dens <- density(cond$median_fc_overall)
  peaks <- which(diff(sign(diff(dens$y))) == -2) + 1
  num_peaks <- length(peaks)
  sd_dist<-c(sd_dist, num_peaks)
  
  
  #sd_dist<-c(sd_dist, sd(cond$median_fc_overall) )
}
p_dist<- ifelse(p_dist<0.05, "S","I")
df_modality<- data.frame( 
  list(conditions=conditions,p_dist=p_dist, d_dist=d_dist, sd_dist=sd_dist )
)

condition_tags
# Load condition tags and clean encodings for plotting labels.
condition_tags<- read_csv("/Data/condition_group.csv", locale = locale(encoding = "UTF-8"))
condition_tags$File_name <- gsub("[^A-Za-z0-9_]", "", condition_tags$File_name)
condition_tags$Clean_Tags <- iconv(condition_tags$Clean_Tags, from = "latin1", to = "UTF-8")

names(df) <- gsub("[^A-Za-z0-9_]", ".", names(df))
df_modality$tag <- condition_tags$Category[match( gsub("[^A-Za-z0-9_]", "", df_modality$conditions)  ,  condition_tags$File_name)]
df_modality$name <- condition_tags$Clean_Tags[match( gsub("[^A-Za-z0-9_]", "", df_modality$conditions)  ,  condition_tags$File_name)]
df_modality<-df_modality[!is.na(df_modality$tag),]

df_modality$tag

df_modality$conditions
unique(df_modality$tag)

# Fixed category order used to keep plots visually consistent.
tag_order <- c(
  "Antibiotics",
  "Antiseptic",
  "Envelope stress",
  "Chemical stress",
  "Metal stress",
  "Environmental stress",
  "Stress combination",
  "Carbon source utilisation",
  "Base"
)

df_modality$name[duplicated(df_modality$name)]

df_modality <- df_modality %>%
  mutate(tag = factor(tag, levels = tag_order)) %>%
  arrange(tag) %>%
  mutate(name = factor(name, levels = name))


# Bar plot of dip statistic for each condition, colored by category.
ggplot(df_modality, aes(x = name, y = d_dist, fill = tag)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = c(
    "Chemical stress" = "#377EB8",            # blue
    "Carbon source utilisation" = "#4DAF4A",  # green
    "Antibiotics" = "#E41A1C",                # red
    "Environmental stress" = "#984EA3",       # purple
    "Antiseptic" = "#FF7F00",                 # orange
    "Stress combination" = "#A65628",         # brown
    "Metal stress" = "#B8860B",               # dark gold
    "Base" = "#999999",                       # grey
    "Envelope stress" = "#17BECF"             # cyan
  )) +
  labs(
    x = "Condition",
    y = "Dip statistic",
    fill = "Condition type"
  ) +
  theme_bw()+
  theme(axis.text.x = element_text( size=5, angle=90, vjust = 0.5,hjust = 1),
        axis.text.y = element_text( size=13, hjust = 1),
        axis.title.x = element_text(color="black", size=15, face="bold"),
        axis.title.y = element_text(color="black", size=15, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        )
  )

# Count per tag and category
# Compute the proportion of unimodal/multimodal conditions per category.
df_prop <- df_modality %>%
  group_by(tag, p_dist) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(tag) %>%
  mutate(prop = n / sum(n))

# Plot
ggplot(df_prop, aes(x = tag, y = prop, fill = p_dist)) +
  geom_col() +
  scale_fill_manual(values = c(
    "S" = "#E41A1C",   # red (multimodal)
    "I" = "#377EB8"    # blue (unimodal)
  )) +
  theme_bw() +
  labs(x = "Condition category",
       y = "Proportion",
       fill = "Modality") +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# -----------------------------------------------------------------------------
# Dip-statistic bar plots with category-level mean reference segments.
# -----------------------------------------------------------------------------
#barplot 
# --- ensure correct ordering ---
df_modality$name <- factor(df_modality$name, levels = unique(df_modality$name))

# --- compute category means + x positions ---
df_mean <- df_modality %>%
  mutate(x = as.numeric(name)) %>%
  group_by(tag) %>%
  summarise(
    mean_d = mean(d_dist, na.rm = TRUE),
    xmin = min(x) - 0.4,
    xmax = max(x) + 0.4,
    .groups = "drop"
  )

# --- colors ---
category_colors <- c(
  "Chemical stress" = "#377EB8",
  "Carbon source utilisation" = "#4DAF4A",
  "Antibiotics" = "#E41A1C",
  "Environmental stress" = "#984EA3",
  "Antiseptic" = "#FF7F00",
  "Stress combination" = "#A65628",
  "Metal stress" = "#B8860B",
  "Base" = "#999999",
  "Envelope stress" = "#17BECF"
)

# --- plot ---
ggplot(df_modality, aes(x = name, y = d_dist, fill = tag)) +
  
  geom_bar(stat = "identity", width = 0.8) +
  
  # 🔥 category-specific horizontal segments
  geom_segment(
    data = df_mean,
    aes(x = xmin, xend = xmax,
        y = mean_d, yend = mean_d,
        color = tag),
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  
  scale_fill_manual(values = category_colors) +
  scale_color_manual(values = category_colors) +
  
  labs(
    x = "Condition",
    y = "Dip statistic",
    fill = "Condition type"
  ) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(color = "black", size = 15, face = "bold"),
    axis.title.y = element_text(color = "black", size = 15, face = "bold")
  )

ggplot(df_modality, aes(x = tag, y = d_dist, fill = tag)) +
  
  # boxplot summarising distribution
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  
  # show individual conditions
  geom_jitter(
    aes(color = tag),
    width = 0.2,
    size = 2,
    alpha = 0.8
  ) +
  
  scale_fill_manual(values = c(
    "Chemical stress" = "#377EB8",
    "Carbon source utilisation" = "#4DAF4A",
    "Antibiotics" = "#E41A1C",
    "Environmental stress" = "#984EA3",
    "Antiseptic" = "#FF7F00",
    "Stress combination" = "#A65628",
    "Metal stress" = "#B8860B",
    "Base" = "#999999",
    "Envelope stress" = "#17BECF"
  )) +
  
  scale_color_manual(values = c(
    "Chemical stress" = "#377EB8",
    "Carbon source utilisation" = "#4DAF4A",
    "Antibiotics" = "#E41A1C",
    "Environmental stress" = "#984EA3",
    "Antiseptic" = "#FF7F00",
    "Stress combination" = "#A65628",
    "Metal stress" = "#B8860B",
    "Base" = "#999999",
    "Envelope stress" = "#17BECF"
  )) +
  
  labs(
    x = "Condition type",
    y = "Dip statistic",
    fill = "Condition type"
  ) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(color="black", size=15, face="bold"),
    axis.title.y = element_text(color="black", size=15, face="bold"),
    legend.position = "none"
  )



df_summary <- df_modality %>%
  group_by(tag, p_dist) %>%
  summarise(n = n(), .groups = "drop")

df_summary <- df_modality %>%
  group_by(tag, p_dist) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(tag) %>%
  mutate(
    total = sum(n),
    prop = n / total
  ) %>%
  ungroup()

df_se <- df_summary %>%
  filter(p_dist == "S") %>%
  mutate(
    se = sqrt(prop * (1 - prop) / total),
    ymin = prop - se,
    ymax = prop + se
  )


ggplot(df_summary, aes(x = tag, y = prop, fill = p_dist)) +
  geom_bar(stat = "identity", width = 0.8) +
  
  geom_errorbar(
    data = df_se,
    aes(x = tag, ymin = ymin, ymax = ymax),
    width = 0.2,
    color = "black"
  ) +
  
  scale_fill_manual(values = c(
    "I" = "#377EB8",
    "S" = "#E41A1C"
  )) +
  
  labs(
    x = "Condition type",
    y = "Relative frequency",
    fill = "Modality"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 11, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(color="black", size=15, face="bold"),
    axis.title.y = element_text(color="black", size=15, face="bold"),
    legend.position = "none"
  )


# -----------------------------------------------------------------------------
# Intraclass correlation coefficient summaries.
# ICC values are extracted from condition summary files and plotted by condition
# and by condition category.
# -----------------------------------------------------------------------------
#ICC
library(diptest)
# Directory containing IRIS colony-measurement output files.
dir_path <- "/Data/assay_results_full"
files <- list.files(dir_path,pattern = "-1_A\\.JPG\\.iris$")
files<-gsub("-1_A.JPG.iris","",files)

conditions<- unique(as.character( sapply(files, function(x) paste0( strsplit(x,"-")[[1]][1], "-",  strsplit(x,"-")[[1]][2]) ) ))  
cond_df<-c()
# Container for one ICC value per condition.
icc_list <-c()

for(i in 1:length(conditions)){
  print(i)
  cond<-read_csv(paste0(dir_path, "/", "condition_ICC_New_",conditions[i], ".csv"),show_col_types = FALSE)
  cond$condition<-conditions[i]
  icc_list <-c(icc_list, cond$ICC_over_plates[1])
}


df_modality<- data.frame( 
  list(conditions=conditions,icc=icc_list )
)



cor(tmp$d_dist[match( df_modality$conditions, tmp$conditions)], df_modality$icc , use ="pairwise.complete.obs")

cor.test(tmp$d_dist[match( df_modality$conditions, tmp$conditions)], df_modality$icc , use ="pairwise.complete.obs")


condition_tags<- read_csv("/Data/condition_group.csv", locale = locale(encoding = "UTF-8"))
condition_tags$File_name <- gsub("[^A-Za-z0-9_]", "", condition_tags$File_name)
condition_tags$Clean_Tags <- iconv(condition_tags$Clean_Tags, from = "latin1", to = "UTF-8")

names(df) <- gsub("[^A-Za-z0-9_]", ".", names(df))
df_modality$tag <- condition_tags$Category[match( gsub("[^A-Za-z0-9_]", "", df_modality$conditions)  ,  condition_tags$File_name)]
df_modality$name <- condition_tags$Clean_Tags[match( gsub("[^A-Za-z0-9_]", "", df_modality$conditions)  ,  condition_tags$File_name)]
df_modality<-df_modality[!is.na(df_modality$tag),]

df_modality$tag

df_modality$conditions
condition_tags$Clean_Tags[which(is.na(condition_tags$Category))]

unique(df_modality$tag)

tag_order <- c(
  "Antibiotics",
  "Antiseptic",
  "Envelope stress",
  "Chemical stress",
  "Metal stress",
  "Environmental stress",
  "Stress combination",
  "Carbon source utilisation",
  "Base"
)

df_modality$name[duplicated(df_modality$name)]

df_modality <- df_modality %>%
  mutate(tag = factor(tag, levels = tag_order)) %>%
  arrange(tag) %>%
  mutate(name = factor(name, levels = name))


# Bar plot of ICC for each condition.
ggplot(df_modality, aes(x = name, y = icc , fill = tag)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = c(
    "Chemical stress" = "#377EB8",            # blue
    "Carbon source utilisation" = "#4DAF4A",  # green
    "Antibiotics" = "#E41A1C",                # red
    "Environmental stress" = "#984EA3",       # purple
    "Antiseptic" = "#FF7F00",                 # orange
    "Stress combination" = "#A65628",         # brown
    "Metal stress" = "#B8860B",               # dark gold
    "Base" = "#999999",                       # grey
    "Envelope stress" = "#17BECF"             # cyan
  )) +
  labs(
    x = "Condition",
    y = "ICC",
    fill = "Condition type"
  ) +
  theme_bw()+
  theme(axis.text.x = element_text( size=5, angle=90, vjust = 0.5,hjust = 1),
        axis.text.y = element_text( size=13, hjust = 1),
        axis.title.x = element_text(color="black", size=15, face="bold"),
        axis.title.y = element_text(color="black", size=15, face="bold"),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold"
        )
  )

# --- ordering ---
df_modality <- df_modality %>%
  mutate(tag = factor(tag, levels = tag_order)) %>%
  arrange(tag) %>%
  mutate(name = factor(name, levels = name))

# --- compute category means + x span ---
df_mean <- df_modality %>%
  mutate(x = as.numeric(name)) %>%
  group_by(tag) %>%
  summarise(
    mean_icc = mean(icc, na.rm = TRUE),
    xmin = min(x) - 0.4,
    xmax = max(x) + 0.4,
    .groups = "drop"
  )

# --- colors ---
category_colors <- c(
  "Chemical stress" = "#377EB8",
  "Carbon source utilisation" = "#4DAF4A",
  "Antibiotics" = "#E41A1C",
  "Environmental stress" = "#984EA3",
  "Antiseptic" = "#FF7F00",
  "Stress combination" = "#A65628",
  "Metal stress" = "#B8860B",
  "Base" = "#999999",
  "Envelope stress" = "#17BECF"
)

# --- plot ---
ggplot(df_modality, aes(x = name, y = icc, fill = tag)) +
  
  geom_bar(stat = "identity", width = 0.8) +
  
  # 🔥 category-specific mean lines
  geom_segment(
    data = df_mean,
    aes(x = xmin, xend = xmax,
        y = mean_icc, yend = mean_icc,
        color = tag),
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  
  scale_fill_manual(values = category_colors) +
  scale_color_manual(values = category_colors) +
  
  labs(
    x = "Condition",
    y = "ICC",
    fill = "Condition type"
  ) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(color = "black", size = 15, face = "bold"),
    axis.title.y = element_text(color = "black", size = 15, face = "bold")
  )






ggplot(df_modality, aes(x = tag, y = icc , fill = tag)) +
  
  # boxplot summarising distribution
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  
  # show individual conditions
  geom_jitter(
    aes(color = tag),
    width = 0.2,
    size = 2,
    alpha = 0.8
  ) +
  
  scale_fill_manual(values = c(
    "Chemical stress" = "#377EB8",
    "Carbon source utilisation" = "#4DAF4A",
    "Antibiotics" = "#E41A1C",
    "Environmental stress" = "#984EA3",
    "Antiseptic" = "#FF7F00",
    "Stress combination" = "#A65628",
    "Metal stress" = "#B8860B",
    "Base" = "#999999",
    "Envelope stress" = "#17BECF"
  )) +
  
  scale_color_manual(values = c(
    "Chemical stress" = "#377EB8",
    "Carbon source utilisation" = "#4DAF4A",
    "Antibiotics" = "#E41A1C",
    "Environmental stress" = "#984EA3",
    "Antiseptic" = "#FF7F00",
    "Stress combination" = "#A65628",
    "Metal stress" = "#B8860B",
    "Base" = "#999999",
    "Envelope stress" = "#17BECF"
  )) +
  
  labs(
    x = "Condition type",
    y = "ICC",
    fill = "Condition type"
  ) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(color="black", size=15, face="bold"),
    axis.title.y = element_text(color="black", size=15, face="bold"),
    legend.position = "none"
  )


# -----------------------------------------------------------------------------
# GWAS gene-summary figures.
# Genes are annotated with AMRFinder classes, linked to condition categories,
# and summarized by odds-ratio direction.
# -----------------------------------------------------------------------------
#GWAS Genes Figure 
# Read pangenome GWAS summary and AMRFinder annotation results.
scatter_gwas <- read_csv("/Data/pangenome_pliotropy.csv",show_col_types = FALSE)
arg_vir<- read_tsv("/Data/AMRFinder_results.tsv",show_col_types = FALSE)

# Map AMR/virulence/stress annotations to GWAS genes.
scatter_gwas$arg_vir <- arg_vir$`Element type`[match(scatter_gwas$Gene, arg_vir$`Contig id` )]
scatter_gwas$arg_vir <- ifelse(grepl("ypothet", scatter_gwas$Annotation),"Hypothetical", scatter_gwas$arg_vir)
scatter_gwas$arg_vir <- ifelse(grepl("hage", scatter_gwas$Annotation),"MGE", scatter_gwas$arg_vir)
scatter_gwas$arg_vir <- ifelse(grepl("ransposon", scatter_gwas$Annotation),"MGE", scatter_gwas$arg_vir)
scatter_gwas$arg_vir <- ifelse(is.na(scatter_gwas$arg_vir),"Other", scatter_gwas$arg_vir)

scatter_gwas_short <- select(scatter_gwas, Gene, Odds_ratio,condition,type,  arg_vir)
scatter_gwas_short$condition <- gsub("results","",scatter_gwas_short$condition)

#conditions<-read_csv("/Data/condition_group.csv",show_col_types = FALSE)
#conditions$File_name <- gsub("[^A-Za-z0-9]+", "", conditions$File_name)

condition_tags<- read_csv("/Data/condition_group.csv", locale = locale(encoding = "UTF-8"))
condition_tags$File_name <- gsub("[^A-Za-z0-9_]", "", condition_tags$File_name)
condition_tags$Clean_Tags <- iconv(condition_tags$Clean_Tags, from = "latin1", to = "UTF-8")

scatter_gwas_short$category<-condition_tags$Category[match(scatter_gwas_short$condition, condition_tags$File_name)]
scatter_gwas_short$Clean_Tags<-condition_tags$Clean_Tags[match(scatter_gwas_short$condition, condition_tags$File_name)]
scatter_gwas_short<-scatter_gwas_short[!is.na(scatter_gwas_short$Clean_Tags),]

#unique(scatter_gwas_short$condition[which(is.na(scatter_gwas_short$Clean_Tags))])

plot_df <- scatter_gwas_short %>%
  mutate(effect = case_when(
    Odds_ratio > 1 ~ "Above 1",
    Odds_ratio < 1 ~ "Below 1",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(effect))



condition_order <- freq_df %>%
  distinct(condition, category) %>%
  arrange(category, condition) %>%
  pull(condition)



# Count GWAS genes by condition, effect direction, gene class, and category.
freq_df <- plot_df %>%
  group_by(condition,Clean_Tags, effect, arg_vir,category) %>%
  summarise(n = n(), .groups = "drop")

freq_df$condition <- factor(freq_df$condition, levels = condition_order)

freq_df$count_new <- 1

freq_df <- freq_df %>%
  mutate(count_mirror =
           ifelse(effect == "Below 1",
                  -count_new,
                  count_new))

unique(freq_df$Clean_Tags)

category_order <- c(
  "Antibiotics",
  "Antiseptic",
  "Envelope stress",
  "Chemical stress",
  "Metal stress",
  "Environmental stress",
  "Stress combination",
  "Carbon source utilisation",
  "Base"
)

freq_df$category <- factor(freq_df$category, levels = category_order)

freq_df <- freq_df |>
  arrange(category, Clean_Tags)

freq_df$Clean_Tags <- factor(freq_df$Clean_Tags,
                             levels = unique(freq_df$Clean_Tags))


ggplot(freq_df,
       aes(x = Clean_Tags,
           y = count_mirror,
           fill = arg_vir)) +
  
  geom_bar(stat = "identity", width = 1) +
  
  geom_hline(yintercept = 0,
             color = "black",
             linewidth = 0.6) +
  
  scale_x_discrete(expand = c(0, 0)) +
  
  scale_fill_manual(values = c(
    AMR = "#8B0000",
    VIRULENCE = "#08306B",
    STRESS = "#006400",
    Other = "grey80",
    Hypothetical = "grey50",
    MGE = "#6A3D9A"
  )) +
  
  facet_grid(~category, scales = "free_x", space = "free_x") +
  
  theme_bw(base_size = 14) +
  
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.5, size = 7),
    axis.ticks.length = unit(2, "pt"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    
    strip.text.x = element_text(angle = 90, face = "bold"),
    strip.background = element_blank(),
    
    panel.spacing.x = unit(0, "pt")
  ) +
  
  labs(
    x = "Condition",
    y = "Gene count (Above ↑  Below ↓)",
    fill = "Gene class"
  )





ggplot(freq_df,
       aes(x = Clean_Tags,
           y = count_mirror,
           fill = arg_vir)) +
  
  geom_col(width = 1) +
  
  geom_hline(yintercept = 0,
             color = "black",
             linewidth = 0.6) +
  
  scale_x_discrete(
    expand = c(0, 0),
    position = "bottom"
  ) +
  
  scale_fill_manual(values = c(
    AMR = "#8B0000",
    VIRULENCE = "#08306B",
    STRESS = "#006400",
    Other = "grey80",
    Hypothetical = "grey50",
    MGE = "#6A3D9A"
  )) +
  
  facet_grid(~category, scales = "free_x", space = "free_x") +
  
  theme_bw(base_size = 14) +
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 7,
      margin = margin(t = 0)
    ),
    
    axis.title.x = element_text(margin = margin(t = 5)),
    
    axis.ticks.x = element_line(),
    axis.ticks.length = unit(2, "pt"),
    
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    
    strip.text.x = element_text(angle = 90, face = "bold"),
    strip.background = element_blank(),
    
    panel.spacing.x = unit(0, "pt")
  ) +
  
  labs(
    x = "Condition",
    y = "Gene count",
    fill = "Gene class"
  )

# -----------------------------------------------------------------------------
# ARG-specific GWAS summary.
# Reads ARG GWAS outputs, filters significant associations, and prepares
# top-gene visualizations.
# -----------------------------------------------------------------------------
#ARGs 
# Folder containing ARG GWAS result CSV files.
folder <- "/Data/ARG"
csv_files <- list.files(path = folder, pattern = "\\.csv$", full.names = TRUE)
output_GWAS_ARG<-c()
# Read each ARG GWAS file and attach the corresponding condition name.
for(i in 1:length(csv_files)){
  trait<-gsub("/Data/ARG/","",csv_files[i])
  trait_nocsv<-gsub(".results.csv","",trait)
  acc_genes<-read_csv(paste0("/Data/ARG/", trait_nocsv, ".results.csv"),show_col_types = FALSE) 
  trait<-gsub(".csv",".",trait)
  trait <- gsub("[^A-Za-z0-9]+", "", trait)
  acc_genes$condition<-trait
  output_GWAS_ARG<-rbind(output_GWAS_ARG, acc_genes)
}

output_GWAS_ARG$type <- ifelse(
  grepl("AMR",output_GWAS_ARG$Gene),"AMR",
  ifelse(grepl("VIRULENCE",output_GWAS_ARG$Gene),"VIRULENCE","OTHERS" )
)

output_GWAS_ARG<- output_GWAS_ARG[ !grepl("rpsJ", output_GWAS_ARG$Gene),]

#write_csv(output_GWAS_pan,"/Data/pangenome_pliotropy.csv" )

output_GWAS_ARG$Gene_short<- as.character( sapply(output_GWAS_ARG$Gene , function(x) strsplit(x,"_")[[1]][1]) )
output_GWAS_ARG_short<- select(output_GWAS_ARG,Gene,Gene_short,Odds_ratio,Benjamini_H_p, Best_pairwise_comp_p, Worst_pairwise_comp_p,condition  )

# Prepare PheWAS-like table with -log10 adjusted p-values and gene classes.
phewas_df <- output_GWAS_ARG_short %>%
  filter(!is.na(Benjamini_H_p)) %>%
  mutate(
    
    logP = -log10(Benjamini_H_p),
    
    gene_class = case_when(
      grepl("_AMR_", Gene) ~ "AMR",
      grepl("_VIRULENCE_", Gene) ~ "VIRULENCE",
      TRUE ~ "Other"
    ),
    
    condition = as.factor(condition)
  )

phewas_df <- phewas_df %>%
  arrange(condition) %>%
  mutate(cond_index = as.numeric(condition))

phewas_df <- phewas_df %>%
  arrange(condition) %>%
  mutate(
    condition=condition,
    cond_index = as.numeric(factor(condition)),
    is_sig = Best_pairwise_comp_p < 0.01
  )

top_genes <- phewas_df %>%
  group_by(condition) %>%
  arrange(desc(logP), .by_group = TRUE) %>%
  distinct(Gene_short, .keep_all = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

# Filter ARG GWAS results using pairwise and FDR-adjusted significance criteria.
output_GWAS_ARG_short_f<- output_GWAS_ARG_short[output_GWAS_ARG_short$Worst_pairwise_comp_p < 0.05,]
output_GWAS_ARG_short_f<- output_GWAS_ARG_short_f[output_GWAS_ARG_short_f$Benjamini_H_p < 0.05,]


top10_per_condition <- output_GWAS_ARG_short_f %>%
  group_by(condition) %>%
  arrange(Benjamini_H_p, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup()

table(top10_per_condition$Gene)


top10_per_condition %>%
  count(condition)

unique(output_GWAS_ARG_short_f$condition)

plot_df$gene_class = case_when(
  grepl("_AMR_", plot_df$Gene) ~ "AMR",
  grepl("_VIRULENCE_", plot_df$Gene) ~ "VIRULENCE",
  TRUE ~ "Other"
)


condition_tags<- read_csv("/Data/condition_group.csv", locale = locale(encoding = "UTF-8"))
condition_tags$File_name <- gsub("[^A-Za-z0-9_]", "", condition_tags$File_name)
condition_tags$Clean_Tags <- iconv(condition_tags$Clean_Tags, from = "latin1", to = "UTF-8")

#scatter_gwas_short$category<-condition_tags$Category[match(scatter_gwas_short$condition, condition_tags$File_name)]
#scatter_gwas_short$Clean_Tags<-condition_tags$Clean_Tags[match(scatter_gwas_short$condition, condition_tags$File_name)]
#scatter_gwas_short<-scatter_gwas_short[!is.na(scatter_gwas_short$Clean_Tags),]

output_GWAS_ARG_short_f$condition <- gsub("results","",output_GWAS_ARG_short_f$condition)
output_GWAS_ARG_short_f$category <- condition_tags$Category[match(output_GWAS_ARG_short_f$condition, condition_tags$File_name)]
output_GWAS_ARG_short_f$Clean_Tags <- condition_tags$Clean_Tags[match(output_GWAS_ARG_short_f$condition, condition_tags$File_name)]
output_GWAS_ARG_short_f <- output_GWAS_ARG_short_f[!is.na(output_GWAS_ARG_short_f$Clean_Tags),]

# prepare top genes and ranking
plot_df <- output_GWAS_ARG_short_f %>%
  group_by(Clean_Tags, condition, category) %>%
  arrange(Benjamini_H_p, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  mutate(
    rank = row_number(),
    OR_direction = ifelse(Odds_ratio > 1, "Above 1", "Below 1")
  ) %>%
  ungroup()

category_order <- c(
  "Antibiotics",
  "Antiseptic",
  "Envelope stress",
  "Chemical stress",
  "Metal stress",
  "Environmental stress",
  "Stress combination",
  "Carbon source utilisation",
  "Base"
)

plot_df <- plot_df[!is.na(plot_df$Clean_Tags),]


plot_df$category <- factor(plot_df$category, levels = category_order)

plot_df$Gene_short <- sub("(_[^_]+){2}$", "", plot_df$Gene)
ggplot(plot_df,
       aes(x = Clean_Tags,
           y = rank,
           fill = log(Odds_ratio+0.000000001))) +
  
  geom_tile(color = "black", linewidth = 0.4) +
  
  geom_text(
    aes(label = Gene_short),
    size = 2,
    fontface = "bold",
    angle = 90,
    vjust = 0.5,
    hjust = 0.5
  ) +
  
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 1,
    name = "Odds ratio"
  ) +
  
  scale_y_reverse(
    breaks = 1:5,
    limits = c(5.5, 0.5),
    expand = c(0,0)
  ) +
  
  scale_x_discrete(expand = c(0,0)) +
  
  facet_grid(~category,
             scales = "free_x",
             space = "free_x") +
  
  labs(
    x = "Condition",
    y = "Top genes rank"
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 7,
      margin = margin(t = 0)
    ),
    
    panel.grid = element_blank(),
    
    strip.text.x = element_text(angle = 90, face = "bold"),
    strip.background = element_blank(),
    
    panel.spacing.x = unit(0, "pt")
  )
plot_df$Odds_ratio


# -----------------------------------------------------------------------------
# Condition-correlation and clustering analysis.
# Pairwise fitness-profile correlations are embedded with t-SNE and clustered.
# -----------------------------------------------------------------------------
#correlation of conditions
Fitness<-read_csv("/Data/absolute_fitness_ML.csv")
colnames(Fitness)[1]<-"Name"


# Initialize the condition-by-condition correlation matrix.
correlation_matrix<-matrix(1, nrow =  dim(Fitness)[2]-1 , ncol =  dim(Fitness)[2]-1)
# Compute pairwise correlations between all condition fitness profiles.
for(i in 2:dim(Fitness)[2]){
  print(i/dim(Fitness)[2])
  for(j in 2:dim(Fitness)[2]){
    correlation_matrix[i-1,j-1]<- cor(Fitness[,i] %>% pull(.), Fitness[,j] %>% pull(.) )
  }
}
colnames(correlation_matrix) <- colnames(Fitness)[-1]
row.names(correlation_matrix) <- colnames(Fitness)[-1]


condition_tags<- read_csv("/Data/condition_group.csv", locale = locale(encoding = "UTF-8"),show_col_types = FALSE)
condition_tags$File_name <- gsub("[^A-Za-z0-9_]", "", condition_tags$File_name)
condition_tags$Clean_Tags <- iconv(condition_tags$Clean_Tags, from = "latin1", to = "UTF-8")

colnames(correlation_matrix)
category <- condition_tags$Category[match(colnames(correlation_matrix), condition_tags$File_name)]
correlation_matrix<- correlation_matrix[which(!is.na(category )),which(!is.na(category ))]
category <- category[complete.cases(category )]

library(Rtsne)
library(ggplot2)
library(dplyr)

# تبدیل همبستگی به distance
# Convert correlations to distances for t-SNE embedding.
dist_matrix <- as.dist(1 - correlation_matrix)

# اجرای t-SNE
set.seed(123)

tsne_out <- Rtsne(
  dist_matrix,
  is_distance = TRUE,
  perplexity = 5,
  verbose = TRUE
)

# ساخت dataframe برای plot
tsne_df <- data.frame(
  x = tsne_out$Y[,1],
  y = tsne_out$Y[,2],
  variable = rownames(correlation_matrix),
  category = category
)

# رنگ‌ها
category_colors <- c(
  "Chemical stress" = "#377EB8",
  "Carbon source utilisation" = "#4DAF4A",
  "Antibiotics" = "#E41A1C",
  "Environmental stress" = "#984EA3",
  "Antiseptic" = "#FF7F00",
  "Stress combination" = "#A65628",
  "Metal stress" = "#B8860B",
  "Base" = "#999999",
  "Envelope stress" = "#17BECF"
)

# رسم



library(mclust)

mclust_res <- Mclust(
  tsne_df[, c("x","y")],
  G = 1:20   # test up to 20 clusters
)
plot(mclust_res, what = "BIC")

# Run model-based clustering
mclust_res <- Mclust(tsne_df[, c("x","y")], G = 1:20)

# Summary
summary(mclust_res)

# Cluster assignments
tsne_df$cluster <- mclust_res$classification

mclust_res$G

ggplot(tsne_df, aes(x, y, color = factor(cluster))) +
  geom_point(size = 4) +
  geom_text(aes(label = variable), vjust = -0.6, size = 3) +
  theme_classic() +
  labs(color = "Mclust cluster")

ggplot(tsne_df, aes(x, y, color = category, label = variable)) +
  geom_point(size = 2, alpha = 0.7, shape = 16) +
  #geom_text(size = 3, vjust = -0.7) +
  stat_ellipse(aes(group = cluster)) +
  scale_color_manual(values = category_colors) +
  theme_classic() +
  labs(
    x = "t-SNE 1",
    y = "t-SNE 2",
    color = "Category"
  )+
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 11, angle = 0),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(color="black", size=15, face="bold"),
    axis.title.y = element_text(color="black", size=15, face="bold"),
    legend.position = "none"
  )

tab <- table(tsne_df$cluster, tsne_df$category)
plot_df <- as.data.frame(tab)
colnames(plot_df) <- c("cluster", "category", "count")
category_colors <- c(
  "Chemical stress" = "#377EB8",
  "Carbon source utilisation" = "#4DAF4A",
  "Antibiotics" = "#E41A1C",
  "Environmental stress" = "#984EA3",
  "Antiseptic" = "#FF7F00",
  "Stress combination" = "#A65628",
  "Metal stress" = "#B8860B",
  "Base" = "#999999",
  "Envelope stress" = "#17BECF"
)

write_csv(tsne_df, "/Data/clusters.csv")
library(ggplot2)

ggplot(plot_df, aes(x = factor(cluster), y = count, fill = category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = category_colors) +
  theme_classic() +
  labs(
    x = "Cluster",
    y = "Number of conditions",
    fill = "Category"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 11, angle = 0),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(color="black", size=15, face="bold"),
    axis.title.y = element_text(color="black", size=15, face="bold"),
    legend.position = "none"
  )

split(tsne_df$variable, tsne_df$cluster)


library(dplyr)

# ensure names align
rownames(correlation_matrix) <- colnames(correlation_matrix)

# map categories
cat_df <- data.frame(
  id = rownames(correlation_matrix),
  category = category
)

# create category-level matrix
categories <- unique(category)

cat_mat <- matrix(NA, nrow = length(categories), ncol = length(categories))
rownames(cat_mat) <- categories
colnames(cat_mat) <- categories

# Average correlations for each pair of condition categories.
for (c1 in categories) {
  for (c2 in categories) {
    
    idx1 <- which(category == c1)
    idx2 <- which(category == c2)
    
    sub_mat <- correlation_matrix[idx1, idx2]
    
    # remove self-correlation bias if same category
    if (c1 == c2) {
      sub_mat <- sub_mat[lower.tri(sub_mat)]
    }
    
    cat_mat[c1, c2] <- mean(sub_mat, na.rm = TRUE)
  }
}
library(corrplot)

# Compute clustering
hc <- hclust(as.dist(1 - cat_mat))

# Plot dendrogram
# Compute clustering ONCE


# clustering
hc <- hclust(as.dist(1 - cat_mat), method = "complete")

# extract order
ord <- hc$order

# 🔥 reverse to match dendrogram visual orientation
ord_rev <- ord

# reorder matrix
cat_mat_ord <- cat_mat[ord_rev, ord_rev]

# plot dendrogram
plot(as.dendrogram(hc), main = "Clustering")

# plot heatmap with SAME order
corrplot(cat_mat_ord,
         method = "color",
         type = "upper",
         order = "original",   # 🔥 do NOT use hclust here
         col = colorRampPalette(c("blue", "white", "red"))(100),
         tl.col = "black",
         tl.cex = 1,
         addCoef.col = "black",
         diag = T)


# -----------------------------------------------------------------------------
# Plasmid/ARG contribution analysis.
# Combines ARG presence/absence summaries with fitness values and fits elastic
# net models to estimate plasmid-associated fitness effects.
# -----------------------------------------------------------------------------
#plasmid cost
library(tidyverse)
plasmid<-read_tsv("/Data/summary_ARG.tsv")
plasmid$`#FILE` <- gsub("abricate_","", plasmid$`#FILE`)
plasmid$`#FILE` <- gsub(".fasta.tsv","", plasmid$`#FILE`)

Fitness<-read_csv("/Data/absolute_fitness_ML.csv")
# Approach 1: using sapply & var (numeric columns)
# Approach: check for >1 unique non-NA value (works generically)
colnames(Fitness)[1]<-"Name"
Fitness<-Fitness[match(plasmid$`#FILE`,Fitness$Name),]

sort(table(Fitness$Name), decreasing =T)

plasmid <- plasmid %>%
  mutate(
    across(
      3:ncol(.),
      ~ ifelse(. == "." | is.na(.), 0, 1)
    )
  )

#plasmid[,3:dim(plasmid)[2]] <- ifelse(plasmid[,3:dim(plasmid)[2]] =="\\.",1,0)
colnames(plasmid)[1]<-"Name"
plasmid<-plasmid[,-2]

plasmid_short<-plasmid[,c(1, match(names(sort(apply(plasmid[,-1], 2, sum), decreasing =T)[1:10]) , colnames(plasmid) ))]
plasmid_short %>%
  summarise(across(-Name, sum, na.rm = TRUE))

# predictors
# Build predictor matrix from selected ARG/plasmid indicators.
X <- plasmid_short %>%
  column_to_rownames("Name") %>%
  as.matrix()



# responses
# Build response matrix from fitness values for matched strains.
Y <- Fitness %>%
  column_to_rownames("Name") %>%
  as.matrix()

# Restrict the plasmid model to a selected sequence type before fitting.
#STs filtering
kleborate<-read_tsv("/Data/tot_results.tsv",show_col_types = FALSE)
kleborate$ST_short <-as.character(sapply(kleborate$ST, function(x) strsplit(x,"-")[[1]][1]))
ST_vector <- kleborate$ST_short[match( row.names(X) , kleborate$strain)]
#Fitness<- Fitness[which(ST_vector=="ST2096"),]
#plasmid <- plasmid[which(ST_vector=="ST2096"),]

X <- X[which(ST_vector=="ST147"),]
Y <- Y [which(ST_vector=="ST147"),]

# ensure same strain order
common <- intersect(rownames(X), rownames(Y))
X <- X[common, ]
Y <- Y[common, ]

library(glmnet)

tot_vals<-c()
# Fit one elastic-net model per condition to estimate plasmid/ARG effects.
for(i in 1:dim(Y)[2]){
  print(i/dim(Y)[2])
  y <- Y[, i]
  
  # fit model
  cvfit <- cv.glmnet(
    X,
    y,
    alpha = 0.5,      # elastic net
    family = "gaussian",
    standardize = TRUE
  )
  
  fit <- glmnet(
    X,
    y,
    alpha = 0.5,
    lambda = cvfit$lambda.min,
    standardize = TRUE
  )
  
  coef_df <- as.data.frame(as.matrix(coef(fit)))
  coef_df <-coef_df[-1,]
  tot_vals<-cbind(tot_vals, coef_df)
}

tot_vals<-data.frame(tot_vals)
colnames(tot_vals)<- colnames(Fitness)[-1]
row.names(tot_vals)<-colnames(plasmid_short)[-1]
tot_vals[1:10,1:10]

dim(tot_vals)


condition_tags<- read_csv("/Data/condition_group.csv", locale = locale(encoding = "UTF-8"),show_col_types = FALSE)
condition_tags$File_name <- gsub("[^A-Za-z0-9_]", "", condition_tags$File_name)
condition_tags$Clean_Tags <- iconv(condition_tags$Clean_Tags, from = "latin1", to = "UTF-8")

colnames(correlation_matrix)


mat <- as.matrix(tot_vals)
colnames(mat) <- condition_tags$Clean_Tags[match(colnames(mat), condition_tags$File_name)]
mat<-mat[,!is.na(colnames(mat))]
category <- condition_tags$Category[match(colnames(mat), condition_tags$Clean_Tags)]


library(pheatmap)

# order matrix by category
# ensure category corresponds to columns
#category <- category[colnames(mat)]

# get ordering index
order_cat <- order(category)
# reorder matrix columns
mat_ordered <- mat[, order_cat]
# reorder category vector
category <- category[order_cat]

# annotation dataframe
annotation <- data.frame(Category = category)
rownames(annotation) <- colnames(mat_ordered)

# category colors
category_colors <- list(
  Category = c(
    "Chemical stress" = "#377EB8",
    "Carbon source utilisation" = "#4DAF4A",
    "Antibiotics" = "#E41A1C",
    "Environmental stress" = "#984EA3",
    "Antiseptic" = "#FF7F00",
    "Stress combination" = "#A65628",
    "Metal stress" = "#B8860B",
    "Base" = "#999999",
    "Envelope stress" = "#17BECF"
  )
)

# force white at zero
rg <- max(abs(mat_ordered))
breaks <- seq(-rg, rg, length.out = 101)

# heatmap
pheatmap(
  mat_ordered,
  color = colorRampPalette(c("blue","white","red"))(100),
  breaks = breaks,
  annotation_col = annotation,
  annotation_colors = category_colors,
  cluster_rows = TRUE,
  cluster_cols = F,
  scale = "none",
  fontsize_row = 8,
  fontsize_col = 6,
  border_color = NA,
  main = "Plasmid contribution to fitness"
)


library(tidyverse)

df <- as.data.frame(mat_ordered)
df$plasmid <- rownames(df)

df_long <- df %>%
  pivot_longer(-plasmid, names_to="condition", values_to="fitness") 
df_long$category <- condition_tags$Category[match( df_long$condition , condition_tags$Clean_Tags)]


ggplot(df_long, aes(plasmid, fitness, fill=category)) +
  geom_boxplot(outlier.size=0.5) +
  scale_fill_manual(values=category_colors$Category) +
  labs(
    x="Condition category",
    y="Fitness effect"
  )+
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(color="black", size=15, face="bold"),
    axis.title.y = element_text(color="black", size=15, face="bold"),
    legend.position = "none"
  )

#write_csv(df_long, "/Data/plasmid_cost_ST147.csv")
# -----------------------------------------------------------------------------
# Compare plasmid-associated fitness effects between ST147 and ST2096.
# -----------------------------------------------------------------------------
#Plasmid count and cost 
tmp_ST147 <- read_csv("/Data/plasmid_cost_ST147.csv")
tmp_ST2096 <- read_csv("/Data/plasmid_cost_ST2096.csv") 

head(tmp_ST147)
head(tmp_ST2096)

library(lme4)
library(lmerTest)

# combine datasets
tmp_ST147$ST <- "ST147"
tmp_ST2096$ST <- "ST2096"

df <- bind_rows(tmp_ST147, tmp_ST2096)

# model
df %>%
  group_by(plasmid) %>%
  group_modify(~ {
    tmp <- .x %>%
      pivot_wider(names_from = ST, values_from = fitness) %>%
      mutate(diff = ST147 - ST2096)
    
    data.frame(
      mean_diff = mean(tmp$diff, na.rm = TRUE),
      p_value = wilcox.test(tmp$diff)$p.value
    )
  })


# -----------------------------------------------------------------------------
# Prediction diagnostic figure: observed versus predicted fitness values for
# the test split, including prediction intervals and the 1:1 reference line.
# -----------------------------------------------------------------------------
#test versus predicted Figure
inp<-read_csv("/Data/prediction_new/new_datamode_predictions_NGBoost_DMSO1_mode4.csv")
inp<-inp[inp$split=="test",]
ggplot(inp, aes(x = actual, y = predicted)) +
  
  # confidence interval vertical lines
  geom_errorbar(
    aes(ymin = lower_95, ymax = upper_95),
    width = 0,
    alpha = 0.4,
    colour = "steelblue"
  ) +
  
  # points
  geom_point(
    aes(colour = split),
    size = 2,
    alpha = 0.4
  ) +
  
  # ideal prediction line
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "red",
    linewidth = 1
  ) +
  
  theme_bw(base_size = 14) +
  
  coord_cartesian(
    xlim = c(0.5, 1.5),
    ylim = c(0.5, 1.5)
  )+
  
  labs(
    title = "",
    x = "Actual fitness",
    y = "Predicted fitness",
    colour = "Dataset split"
  )



# -----------------------------------------------------------------------------
# Machine-learning performance summaries across models, feature modes, and traits.
# -----------------------------------------------------------------------------
#Machine Learning plots
# List all model-result files for the prediction benchmark.
all_files <- list.files(path = "/Data/results_prediction",
                        full.names = TRUE,
                        recursive = F)


all_files_short <- list.files(path = "/Data/results_prediction",
                              full.names = F,
                              recursive = F)

all_files_short <- gsub("new_results_","",all_files_short)


unitigs_files<-list.files(path = "/Data/prediction_new",
                          full.names = F,
                          pattern = "new_ML_unitig*",
                          recursive = F)

unitigs_files<-gsub("new_ML_unitig_input_","",unitigs_files)
unitigs_files<-gsub(".csv","",unitigs_files)

snp_files<-list.files(path = "/Data/prediction_new",
                      full.names = F,
                      pattern = "new_ML_SNP_input_*",
                      recursive = F)

snp_files<-gsub("new_ML_SNP_input_","",snp_files)
snp_files<-gsub(".csv","",snp_files)



model<- as.character(sapply( all_files_short , function(x) strsplit(x,"_")[[1]][4]))
condition<- as.character(sapply( all_files_short , function(x) strsplit(x,"_")[[1]][5]))
mode <- as.character(sapply( all_files_short , function(x) strsplit(x,"_")[[1]][6]))
mode <- gsub(".csv","", mode)

# Remove results for feature modes where required input files are unavailable.
ex_list_unitig <- which(!condition %in% unitigs_files & mode == "mode4")
ex_list_snps <- which(!condition %in% snp_files & mode == "mode3")

all_files_short<-all_files_short[-c(ex_list_unitig, ex_list_snps)]
all_files <-all_files[-c(ex_list_unitig, ex_list_snps)]

model<- as.character(sapply( all_files_short , function(x) strsplit(x,"_")[[1]][4]))
condition<- as.character(sapply( all_files_short , function(x) strsplit(x,"_")[[1]][5]))
mode <- as.character(sapply( all_files_short , function(x) strsplit(x,"_")[[1]][6]))
mode <- gsub(".csv","", mode)


# Containers for held-out performance metrics from each result file.
rho_test<-c()
cov95_test<-c()
width_test<-c()

for(i in 1:length(all_files )){
  print(i/length(all_files ))
  tmp<- read_csv(all_files[i],show_col_types = FALSE)
  tmp<- tmp[tmp$split=="test",]
  rho_test<-c(rho_test,tmp$rho )
  cov95_test<-c(cov95_test, tmp$cov95 )
  width_test<-c(width_test, tmp$width )
}

performance_df<- data.frame(
  model=model,
  condition=condition,
  mode=mode,
  rho_test=rho_test,
  cov95_test=cov95_test,
  width_test=width_test
)


df <- performance_df

# ensure mode ordered correctly
df <- df %>%
  mutate(mode = factor(mode, levels = c("mode1","mode2","mode3","mode4","mode5")))

# ensure factors
df_results <- df %>%
  mutate(
    mode  = factor(mode),
    model = factor(model)
  )

# plot 1 
ggplot(df_results, aes(x = model , y = rho_test, fill = mode)) +
  # violin distribution
  geom_violin(
    position = position_dodge(width = 0.8),
    alpha = 0.5,
    trim = FALSE
  ) +
  # boxplot summary
  geom_boxplot(
    width = 0.2,
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    alpha = 0.08
  ) +
  # individual points (each condition)
  geom_jitter(
    aes(color = mode),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    size = 1,
    alpha = 0.06
  ) +
  
  labs(
    x = "Feature Mode",
    y = "Spearman Correlation (ρ)",
    fill = "Model",
    color = "Model"
  ) +
 theme_bw() +
  theme(
    axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(color="black", size=15, face="bold"),
    axis.title.y = element_text(color="black", size=15, face="bold"),
    legend.position = "none"
  )


# plot 2 coverage 
ggplot(df_results, aes(x = model , y = cov95_test, fill = mode)) +
  # violin distribution
  geom_violin(
    position = position_dodge(width = 0.8),
    alpha = 0.5,
    trim = FALSE
  ) +
  # boxplot summary
  geom_boxplot(
    width = 0.2,
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    alpha = 0.08
  ) +
  # individual points (each condition)
  geom_jitter(
    aes(color = mode),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    size = 1,
    alpha = 0.06
  ) +
  
  labs(
    x = "Feature Mode",
    y = "Coverage 95%",
    fill = "Model",
    color = "Model"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(color="black", size=15, face="bold"),
    axis.title.y = element_text(color="black", size=15, face="bold"),
    legend.position = "none"
  )


# Directory containing IRIS colony-measurement output files.
dir_path <- "/Data/assay_results_full"
files <- list.files(dir_path,pattern = "condition_ICC_New_")
ICC_vals<-c()
sd_vals<-c()

condition<-c()
for(i in 1:length(files)){
  tmp<-read_csv(paste0(dir_path, "/", files[i]),show_col_types = FALSE)
  tmp_tag<- gsub("condition_ICC_New_","",files[i])
  tmp_tag<- gsub(".csv","",tmp_tag)
  
  condition<-c(condition, gsub("[^A-Za-z0-9]", "", tmp_tag))
  ICC_vals<-c(ICC_vals, tmp$ICC_over_plates[1] )
  sd_vals<-c(sd_vals,sd(tmp$median_fc_overall))
}
ICC_df=data.frame(
  condition=condition,
  ICC_vals=ICC_vals,
  sd_vals=sd_vals
)

performance_df$ICC_val<- ICC_df$ICC_vals[match(performance_df$condition, ICC_df$condition)]
performance_df$rho_norm <- performance_df$rho_test / sqrt(performance_df$ICC_val)
hist(performance_df$rho_norm)


df_results$sd <- ICC_df$sd_vals[match( df_results$condition,  ICC_df$condition)]
df_results$width_ratio <- df_results$width_test/df_results$sd 
# plot 3 
ggplot(df_results, aes(x = model , y = width_ratio , fill = mode)) +
  # violin distribution
  geom_violin(
    position = position_dodge(width = 0.8),
    alpha = 0.5,
    trim = FALSE
  ) +
  # boxplot summary
  geom_boxplot(
    width = 0.2,
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    alpha = 0.08
  ) +
  # individual points (each condition)
  geom_jitter(
    aes(color = mode),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    size = 1,
    alpha = 0.06
  ) +
  
  labs(
    x = "Feature Mode",
    y = "Width normalized by sd",
    fill = "Model",
    color = "Model"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(color="black", size=15, face="bold"),
    axis.title.y = element_text(color="black", size=15, face="bold"),
    legend.position = "none"
  )


# -----------------------------------------------------------------------------
# Overfitting diagnostics: compare held-out performance against validation and
# training summaries using performance ratios.
# -----------------------------------------------------------------------------
#overfitting 
rho_train_test<-c()
rho_val_test<-c()

for(i in 1:length(all_files )){
  print(i/length(all_files ))
  tmp_test<- read_csv(all_files[i],show_col_types = FALSE)
  tmp_test<- tmp_test[tmp_test$split=="test",]
  
  tmp_val<- read_csv(all_files[i],show_col_types = FALSE)
  tmp_val<- tmp_val[tmp_val$split=="val_cv_mean",]
  
  tmp_train<- read_csv(all_files[i],show_col_types = FALSE)
  tmp_train<- tmp_train[tmp_train$split=="train_cv_mean",]
  
  rho_train_test<-c(rho_train_test,tmp_test$rho/tmp_val$rho )
  rho_val_test<-c(rho_val_test, tmp_test$rho/tmp_train$rho )
}



length(condition)
length(rho_train_test)

model<- as.character(sapply( all_files_short  , function(x) strsplit(x,"_")[[1]][4]))
condition<- as.character(sapply( all_files_short   , function(x) strsplit(x,"_")[[1]][5]))
mode <- as.character(sapply( all_files_short   , function(x) strsplit(x,"_")[[1]][6]))
mode <- gsub(".csv","", mode)


performance_df<- data.frame(
  model=model,
  condition=condition,
  mode=mode,
  rho_train_test=rho_train_test,
  rho_val_test=rho_val_test
)
performance_df<-performance_df[performance_df$rho_train_test>0,]
which(performance_df$rho_train_test<0)
ggplot(performance_df, aes(x = model , y = rho_train_test , fill = mode)) +
  # violin distribution
  geom_violin(
    position = position_dodge(width = 0.8),
    alpha = 0.5,
    trim = FALSE
  ) +
  # boxplot summary
  geom_boxplot(
    width = 0.2,
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    alpha = 0.08
  ) +
  ylim(range(0,4))+
  # individual points (each condition)
  geom_jitter(
    aes(color = mode),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    size = 1,
    alpha = 0.06
  ) +
  
  labs(
    x = "Feature Mode",
    y = "Test over Train",
    fill = "Model",
    color = "Model"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(color="black", size=15, face="bold"),
    axis.title.y = element_text(color="black", size=15, face="bold"),
    legend.position = "none"
  )


#val test

performance_df<- data.frame(
  model=model,
  condition=condition,
  mode=mode,
  rho_train_test=rho_train_test,
  rho_val_test=rho_val_test
)
performance_df<-performance_df[performance_df$rho_val_test>0,]
which(performance_df$rho_val_test<0)
ggplot(performance_df, aes(x = model , y = rho_val_test , fill = mode)) +
  # violin distribution
  geom_violin(
    position = position_dodge(width = 0.8),
    alpha = 0.5,
    trim = FALSE
  ) +
  # boxplot summary
  geom_boxplot(
    width = 0.2,
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    alpha = 0.08
  ) +
  ylim(range(0,2))+
  # individual points (each condition)
  geom_jitter(
    aes(color = mode),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    size = 1,
    alpha = 0.06
  ) +
  
  labs(
    x = "Feature Mode",
    y = "Test over Validation",
    fill = "Model",
    color = "Model"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(color="black", size=15, face="bold"),
    axis.title.y = element_text(color="black", size=15, face="bold"),
    legend.position = "none"
  )





# -----------------------------------------------------------------------------
# Best-model-per-condition summaries.
# Select the highest-performing model per condition and visualize normalized
# performance, coverage, interval width, and model preferences.
# -----------------------------------------------------------------------------
#new plotting prediction 

all_files <- list.files(path = "/Data/results_prediction",
                        full.names = TRUE,
                        recursive = F)


all_files_short <- list.files(path = "/Data/results_prediction",
                              full.names = F,
                              recursive = F)

all_files_short <- gsub("new_results_","",all_files_short)


unitigs_files<-list.files(path = "/Data/prediction_new",
                          full.names = F,
                          pattern = "new_ML_unitig*",
                          recursive = F)

unitigs_files<-gsub("new_ML_unitig_input_","",unitigs_files)
unitigs_files<-gsub(".csv","",unitigs_files)

snp_files<-list.files(path = "/Data/prediction_new",
                      full.names = F,
                      pattern = "new_ML_SNP_input_*",
                      recursive = F)

snp_files<-gsub("new_ML_SNP_input_","",snp_files)
snp_files<-gsub(".csv","",snp_files)

model<- as.character(sapply( all_files_short , function(x) strsplit(x,"_")[[1]][4]))
condition<- as.character(sapply( all_files_short , function(x) strsplit(x,"_")[[1]][5]))
mode <- as.character(sapply( all_files_short , function(x) strsplit(x,"_")[[1]][6]))
mode <- gsub(".csv","", mode)

ex_list_unitig <- which(!condition %in% unitigs_files & mode == "mode4")
ex_list_snps <- which(!condition %in% snp_files & mode == "mode3")

all_files_short<-all_files_short[-c(ex_list_unitig, ex_list_snps)]
all_files <-all_files[-c(ex_list_unitig, ex_list_snps)]


model<- as.character(sapply( all_files_short , function(x) strsplit(x,"_")[[1]][4]))
condition<- as.character(sapply( all_files_short , function(x) strsplit(x,"_")[[1]][5]))
mode <- as.character(sapply( all_files_short , function(x) strsplit(x,"_")[[1]][6]))
mode <- gsub(".csv","", mode)


rho_test<-c()
cov95_test<-c()
width_test<-c()

for(i in 1:length(all_files )){
  print(i/length(all_files ))
  tmp<- read_csv(all_files[i],show_col_types = FALSE)
  tmp<- tmp[tmp$split=="test",]
  rho_test<-c(rho_test,tmp$rho )
  cov95_test<-c(cov95_test, tmp$cov95 )
  width_test<-c(width_test, tmp$width )
}

performance_df<- data.frame(
  model=model,
  condition=condition,
  mode=mode,
  rho_test=rho_test,
  cov95_test=cov95_test,
  width_test=width_test
)

# Select the model/mode with maximum test correlation for each condition.
best_per_condition <- performance_df %>%
  group_by(condition) %>%
  slice_max(order_by = rho_test, n = 1, with_ties = FALSE) %>%
  ungroup()

table(best_per_condition$model)
table(best_per_condition$mode)


boxplot(performance_df$rho_test ~ performance_df$model)
boxplot(performance_df$cov95_test ~ performance_df$model)


# Directory containing IRIS colony-measurement output files.
dir_path <- "/Data/assay_results_full"
files <- list.files(dir_path,pattern = "condition_ICC_New_")
ICC_vals<-c()
sd_vals<-c()

condition<-c()
for(i in 1:length(files)){
  tmp<-read_csv(paste0(dir_path, "/", files[i]),show_col_types = FALSE)
  tmp_tag<- gsub("condition_ICC_New_","",files[i])
  tmp_tag<- gsub(".csv","",tmp_tag)
  
  condition<-c(condition, gsub("[^A-Za-z0-9]", "", tmp_tag))
  ICC_vals<-c(ICC_vals, tmp$ICC_over_plates[1] )
  sd_vals<-c(sd_vals,sd(tmp$median_fc_overall))
}
ICC_df=data.frame(
  condition=condition,
  ICC_vals=ICC_vals,
  sd_vals=sd_vals
)

best_per_condition$ICC_val<- ICC_df$ICC_vals[match(best_per_condition$condition, ICC_df$condition)]
best_per_condition$rho_norm <- best_per_condition$rho_test / sqrt(best_per_condition$ICC_val)
best_per_condition$sd_vals <- ICC_df$sd_vals[match(best_per_condition$condition, ICC_df$condition)]
best_per_condition$width_test_normalized <- best_per_condition$width_test / best_per_condition$sd_vals

hist(best_per_condition$rho_norm)


condition_tags<- read_csv("/Data/condition_group.csv", locale = locale(encoding = "UTF-8"),show_col_types = FALSE)
condition_tags$File_name <- gsub("[^A-Za-z0-9_]", "", condition_tags$File_name)
condition_tags$Clean_Tags <- iconv(condition_tags$Clean_Tags, from = "latin1", to = "UTF-8")

colnames(mat)  <- condition_tags$Clean_Tags[match(colnames(mat), condition_tags$File_name)]

best_per_condition$Clean_Tags <- condition_tags$Clean_Tags[match( best_per_condition$condition, condition_tags$File_name)]
best_per_condition$Category <- condition_tags$Category[match( best_per_condition$condition, condition_tags$File_name)]

best_per_condition <- best_per_condition[!is.na(best_per_condition$Clean_Tags),]


category_colors <- list(
  Category = c(
    "Chemical stress" = "#377EB8",
    "Carbon source utilisation" = "#4DAF4A",
    "Antibiotics" = "#E41A1C",
    "Environmental stress" = "#984EA3",
    "Antiseptic" = "#FF7F00",
    "Stress combination" = "#A65628",
    "Metal stress" = "#B8860B",
    "Base" = "#999999",
    "Envelope stress" = "#17BECF"
  )
)

ggplot(best_per_condition,
       aes(x = Category,
           y = rho_norm,
           fill = Category)) +
  
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  
  geom_jitter(aes(color = Category),
              width = 0.2,
              alpha = 0.6,
              size = 2) +
  
  scale_fill_manual(values = category_colors$Category) +
  scale_color_manual(values = category_colors$Category) +
  
  theme_bw(base_size = 14) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  
  labs(
    x = "Condition category",
    y = expression(rho[test])
  )


ggplot(best_per_condition,
       aes(x = Category,
           y = cov95_test,
           fill = Category)) +
  
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  
  geom_jitter(aes(color = Category),
              width = 0.2,
              alpha = 0.6,
              size = 2) +
  
  scale_fill_manual(values = category_colors$Category) +
  scale_color_manual(values = category_colors$Category) +
  
  theme_bw(base_size = 14) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  
  labs(
    x = "Condition category",
    y = "Coverage 95%"
  )


ggplot(best_per_condition,
       aes(x = Category,
           y = width_test_normalized,
           fill = Category)) +
  
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  
  geom_jitter(aes(color = Category),
              width = 0.2,
              alpha = 0.6,
              size = 2) +
  
  scale_fill_manual(values = category_colors$Category) +
  scale_color_manual(values = category_colors$Category) +
  
  theme_bw(base_size = 14) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  
  labs(
    x = "Condition category",
    y = "Prediction Width"
  )


best_per_condition
table(best_per_condition$Category,best_per_condition$model)


plot_df <- best_per_condition %>%
  count(Category, model)


model_colors <- c(
  "NGBoost" = "#0072B2",  # blue
  "QBoost" = "#E69F00",  # orange
  "QRF" = "#009E73",  # bluish green
  "QRegression" = "#F55E00"   # vermillion
)

ggplot(plot_df,
       aes(x = Category,
           y = n,
           fill = model)) +
  
  geom_bar(stat = "identity") +
  
  scale_fill_manual(values = model_colors) +
  
  theme_bw(base_size = 14) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  
  labs(
    x = "Condition category",
    y = "Number of conditions",
    fill = "Model"
  )

ggplot(plot_df,
       aes(x = Category,
           y = n,
           fill = model)) +
  
  geom_bar(stat = "identity", position = "fill") +
  
  scale_fill_manual(values = model_colors) +
  
  theme_bw(base_size = 14) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  
  labs(
    x = "Condition category",
    y = "Number of conditions",
    fill = "Model"
  )

# -----------------------------------------------------------------------------
# Condition-level best-model bar plots with category mean performance overlays.
# -----------------------------------------------------------------------------
#Last ML plots

category_colors <- c(
  "Chemical stress" = "#377EB8",
  "Carbon source utilisation" = "#4DAF4A",
  "Antibiotics" = "#E41A1C",
  "Environmental stress" = "#984EA3",
  "Antiseptic" = "#FF7F00",
  "Stress combination" = "#A65628",
  "Metal stress" = "#B8860B",
  "Base" = "#999999",
  "Envelope stress" = "#17BECF"
)

plot_df <- best_per_condition %>%
  arrange(Category, Clean_Tags)

plot_df$Clean_Tags <- factor(plot_df$Clean_Tags,
                             levels = unique(plot_df$Clean_Tags))


category_colors <- c(
  "Chemical stress" = "#377EB8",
  "Carbon source utilisation" = "#4DAF4A",
  "Antibiotics" = "#E41A1C",
  "Environmental stress" = "#984EA3",
  "Antiseptic" = "#FF7F00",
  "Stress combination" = "#A65628",
  "Metal stress" = "#B8860B",
  "Base" = "#999999",
  "Envelope stress" = "#17BECF"
)

plot_df <- best_per_condition %>%
  arrange(Category, Clean_Tags)

plot_df$Clean_Tags <- factor(plot_df$Clean_Tags,
                             levels = unique(plot_df$Clean_Tags))

# compute mean per category
mean_df <- plot_df %>%
  group_by(Category) %>%
  summarise(mean_rho = mean(rho_norm, na.rm = TRUE))

# get x positions
plot_df$Clean_Tags <- factor(plot_df$Clean_Tags,
                             levels = unique(plot_df$Clean_Tags))

# category ranges on x-axis
ranges <- plot_df %>%
  mutate(x = as.numeric(Clean_Tags)) %>%
  group_by(Category) %>%
  summarise(xmin = min(x) - 0.4,
            xmax = max(x) + 0.4)

mean_df <- left_join(mean_df, ranges, by = "Category")

ggplot(plot_df,
       aes(x = Clean_Tags,
           y = rho_norm,
           fill = Category)) +
  
  geom_col(width = 0.8) +
  
  scale_fill_manual(values = category_colors) +
  
  geom_segment(data = mean_df,
               aes(x = xmin,
                   xend = xmax,
                   y = mean_rho,
                   yend = mean_rho,
                   color = Category),
               linewidth = 0.5,
               inherit.aes = FALSE) +
  
  scale_color_manual(values = category_colors, guide = "none") +
  
  theme_bw(base_size = 14) +
  
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5.5),
    panel.grid.minor = element_blank()
  ) +
  
  labs(
    x = "Condition",
    y = expression(rho[norm]),
    fill = "Category"
  )

#Tree show an example 



# -----------------------------------------------------------------------------
# Training-ratio sensitivity analysis.
# Summarize model performance as a function of training-data ratio.
# -----------------------------------------------------------------------------
#Ratio plot

#ratio plotting 
files <- list.files(path = "/Data/prediction_new",
                    pattern = "new_results_ratio_*",
                    full.names = T) 
results_tot<-c()
for(j in 1:length(files)){
  print(j/length(files))
  tmp<-read_csv(files[j],show_col_types = FALSE)
  results_tot<-rbind(results_tot,tmp )
}

# --------------------------------------------------
# 1) Summarize by ratio & split
# --------------------------------------------------
# Summarize mean and standard deviation of metrics by ratio and split.
results_summary <- results_tot %>%
  group_by(ratio, split) %>%
  summarize(
    mean_rho  = mean(rho),
    sd_rho    = sd(rho),
    mean_cov95 = mean(cov95),
    sd_cov95   = sd(cov95),
    mean_width = mean(width),
    sd_width   = sd(width),
    .groups = "drop"
  )

# View summary
print(results_summary)

# --------------------------------------------------
# 2) Plot mean across targets with ribbon
# --------------------------------------------------

# Plot for rho
ggplot(results_summary, aes(x = ratio*5, y = mean_rho, color = split, fill = split)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean_rho - sd_rho, ymax = mean_rho + sd_rho),
              alpha = 0.2, linetype = "dashed") +
  ylim(0,1) +
  scale_x_continuous(
    limits = c(5, 100),
    breaks = seq(5, 100, by = 5)
  )+
  labs(
    x = "Ratio (%)",
    y = expression("Mean Spearman " * rho)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(color="black", size=15, face="bold"),
    axis.title.y = element_text(color="black", size=15, face="bold"),
    legend.position = "none"
  )


# -----------------------------------------------------------------------------
# Stratified versus non-stratified split comparison.
# -----------------------------------------------------------------------------
#barplots stratified non-stratified 

#stratified 
files <- list.files(path = "/Data/prediction_new",
                    pattern = "new_stratified_results*",
                    full.names = T) 

results_tot_stratified<-c()
for(j in 1:length(files)){
  print(j/length(files))
  tmp<-read_csv(files[j],show_col_types = FALSE)
  results_tot_stratified<-rbind(results_tot_stratified,tmp )
}

results_tot_stratified<-results_tot_stratified[results_tot_stratified$split=="test",]
results_tot_stratified$run<-"stratified"

ggplot(results_tot_stratified, aes(x = target, y = rho)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  labs(
    title = "Boxplots of rho by Target",
    x = "Target",
    y = "rho"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )



files <- list.files(path = "/Data/prediction_new",
                    pattern = "nonstratified_results*",
                    full.names = T) 

results_tot_non_stratified<-c()
for(j in 1:length(files)){
  print(j/length(files))
  tmp<-read_csv(files[j],show_col_types = FALSE)
  results_tot_non_stratified<-rbind(results_tot_non_stratified,tmp )
}

results_tot_non_stratified<-results_tot_non_stratified[results_tot_non_stratified$split=="test",]
results_tot_non_stratified$run<-"nonstratified"

ggplot(results_tot_non_stratified, aes(x = target, y = rho)) +
  geom_boxplot(fill = "red", color = "red4") +
  labs(
    title = "Boxplots of rho by Target",
    x = "Target",
    y = "rho"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


results_tot <- rbind(results_tot_stratified, results_tot_non_stratified)


condition_tags<- read_csv("/Data/condition_group.csv", locale = locale(encoding = "UTF-8"),show_col_types = FALSE)
condition_tags$File_name <- gsub("[^A-Za-z0-9_]", "", condition_tags$File_name)
condition_tags$Clean_Tags <- iconv(condition_tags$Clean_Tags, from = "latin1", to = "UTF-8")

results_tot$target <- gsub("[^A-Za-z0-9_]", "", results_tot$target)
results_tot$Clean_Tags <- condition_tags$Clean_Tags[match(results_tot$target, condition_tags$File_name)]
results_tot$Category <- condition_tags$Category[match(results_tot$target, condition_tags$File_name)]

results_tot <- results_tot[!is.na(results_tot$Clean_Tags),]

head(results_tot)


# Order Clean_Tags by Category
results_tot <- results_tot %>%
  arrange(Category, Clean_Tags) %>%
  mutate(Clean_Tags = factor(Clean_Tags, levels = unique(Clean_Tags)))

ggplot(results_tot, aes(x = Clean_Tags, y = rho, fill = run)) +
  geom_boxplot(position = position_identity(), alpha = 0.32) +
  scale_fill_manual(values = c("stratified" = "lightblue",
                               "nonstratified" = "salmon")) +
  labs(
    x = "Condition",
    y = expression(rho),
    fill = "Run"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(color="black", size=15, face="bold"),
    axis.title.y = element_text(color="black", size=15, face="bold"),
    legend.position = "none"
  )


# compute median rho per condition and run
median_df <- results_tot %>%
  group_by(Clean_Tags, Category, run) %>%
  summarise(median_rho = median(rho), .groups = "drop")

# reshape to compute difference
# Compute median performance difference between non-stratified and stratified runs.
diff_df <- median_df %>%
  pivot_wider(names_from = run, values_from = median_rho) %>%
  mutate(diff =  nonstratified - stratified)

# keep categories grouped
diff_df <- diff_df %>%
  arrange(Category, Clean_Tags) %>%
  mutate(Clean_Tags = factor(Clean_Tags, levels = unique(Clean_Tags)))

# plot

category_colors <- list(
  Category = c(
    "Chemical stress" = "#377EB8",
    "Carbon source utilisation" = "#4DAF4A",
    "Antibiotics" = "#E41A1C",
    "Environmental stress" = "#984EA3",
    "Antiseptic" = "#FF7F00",
    "Stress combination" = "#A65628",
    "Metal stress" = "#B8860B",
    "Base" = "#999999",
    "Envelope stress" = "#17BECF"
  )
)

diff_df$Category <- factor(
  diff_df$Category,
  levels = names(category_colors$Category)
)

ggplot(diff_df, aes(x = Clean_Tags, y = diff, fill = Category)) +
  geom_col() +
  scale_fill_manual(values = category_colors$Category) +
  labs(
    title = "Difference in Median ρ (Non-stratified − Stratified)",
    x = "Condition",
    y = expression(Delta~rho)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 9, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(color = "black", size = 15, face = "bold"),
    axis.title.y = element_text(color = "black", size = 15, face = "bold"),
    legend.position = "none"
  )


library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

# --- min–max scaling with sqrt stretch for visibility
# Min-max scale each condition and apply square-root stretch for visibility.
mat_norm <- apply(mat_ordered, 2, function(x) {
  x_scaled <- (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  sqrt(x_scaled)  # expands lower range so differences are more visible
})

mat_norm <- as.matrix(mat_norm)
rownames(mat_norm) <- rownames(mat_ordered)

# category aligned with columns
Category <- condition_tags$Category[
  match(colnames(mat_norm), condition_tags$Clean_Tags)
]





# convert to long format
df_long <- as.data.frame(mat_norm) %>%
  rownames_to_column("ST") %>%
  pivot_longer(-ST, names_to = "condition", values_to = "fitness")

# attach category
df_long$Category <- Category[
  match(df_long$condition, colnames(mat_norm))
]

# mean per ST per category
cat_mean <- df_long %>%
  group_by(ST, Category) %>%
  summarise(mean_fitness = mean(fitness, na.rm = TRUE), .groups = "drop")

# preserve ST order
cat_mean$ST <- factor(cat_mean$ST, levels = ordered_rows)

# preserve category color order
cat_mean$Category <- factor(
  cat_mean$Category,
  levels = names(category_colors$Category)
)

library(dplyr)
library(ggplot2)

# mean and sd per ST per category
cat_mean <- df_long %>%
  group_by(ST, Category) %>%
  summarise(
    mean_fitness = mean(fitness, na.rm = TRUE),
    sd_fitness   = sd(fitness, na.rm = TRUE),
    .groups = "drop"
  )

# preserve ST order
cat_mean$ST <- factor(cat_mean$ST, levels = ordered_rows)

# preserve category order
cat_mean$Category <- factor(cat_mean$Category,
                            levels = names(category_colors$Category))


ggplot(cat_mean,
       aes(x = ST,
           y = mean_fitness,
           fill = Category)) +
  
  geom_col(width = 0.75) +
  
  geom_errorbar(aes(ymin =  pmax(mean_fitness - sd_fitness, 0),
                    ymax = mean_fitness + sd_fitness),
                width = 0.25,
                linewidth = 0.5) +
  
  scale_fill_manual(values = category_colors$Category) +
  
  facet_wrap(~Category, nrow = 1, scales = "free_x") +
  
  coord_flip() +
  
  labs(
    x = "ST",
    y = "Scaled Mean Fitness"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(color = "black", size = 15, face = "bold"),
    axis.title.y = element_text(color = "black", size = 15, face = "bold"),
    legend.position = "none"
  )
#Details machine learning performance 


# -----------------------------------------------------------------------------
# Mantel test between genomic distance and fitness-profile distance.
# -----------------------------------------------------------------------------
#Mantel test 
kmer<-read_csv("/Data/total_distance_matrix.csv")
fitness<-read_csv("/Data/absolute_fitness_ML.csv")

shared<- intersect(colnames(kmer),fitness$Name )
kmer<-kmer[match(shared,colnames(kmer)),match(shared,colnames(kmer))]

fitness<-fitness[match(shared, fitness$Name),]
corr_mat<-matrix(0, nrow = length(shared), ncol = length(shared))

# Build pairwise strain fitness-correlation matrix.
for(i in 1:length(shared)){
  print(i/length(shared))
  for(j in 1:length(shared)){
    corr_mat[i,j] <- cor(as.numeric(fitness[i, -1]),as.numeric(fitness[j, -1]))
  }
}

library(vegan)
kmer_mat <- as.matrix(kmer)
dist_kmer <- dist(kmer_mat, method = "euclidean")
dist_corr <- as.dist(1 - corr_mat)
mantel_res <- mantel(xdis = dist_kmer, ydis = dist_corr,
                     method = "pearson", permutations = 99)

mantel_res$signif

# -----------------------------------------------------------------------------
# Virulence/resistance association models.
# Fit per-condition linear models using virulence score, resistance score, and
# their interaction as predictors of fitness.
# -----------------------------------------------------------------------------
#intercation virulence rsistance
library(pheatmap)
fitness_mat<-read_csv("/Data/absolute_fitness_ML.csv")
# Approach 1: using sapply & var (numeric columns)
# Approach: check for >1 unique non-NA value (works generically)
colnames(fitness_mat)[1]<-"Name"


kleborate<-read_tsv("/Data/tot_results.tsv",show_col_types = FALSE)
kleborate$ST_short <-as.character(sapply(kleborate$ST, function(x) strsplit(x,"-")[[1]][1]))
kleborate$strain <- toupper(kleborate$strain)

shared<- intersect(fitness_mat$Name, kleborate$strain)
fitness_mat<-fitness_mat[match(shared, fitness_mat$Name), ]
kleborate<-kleborate[ match(shared, kleborate$strain),]


condition_tags<- read_csv("/Data/condition_group.csv", locale = locale(encoding = "UTF-8"),show_col_types = FALSE)
condition_tags$File_name <- gsub("[^A-Za-z0-9_]", "", condition_tags$File_name)
condition_tags$Clean_Tags <- iconv(condition_tags$Clean_Tags, from = "latin1", to = "UTF-8")


colnames(fitness_mat)  <- condition_tags$Clean_Tags[match(colnames(fitness_mat), condition_tags$File_name)]
fitness_mat<- fitness_mat[,!is.na(colnames(fitness_mat))]

Category <- condition_tags$Category[match(colnames(fitness_mat), condition_tags$Clean_Tags)]


output<-c()
for(i in 2:dim(fitness_mat)[2]){
  model<-lm( fitness_mat[,i] %>% pull(.) ~  kleborate$virulence_score + kleborate$resistance_score + kleborate$virulence_score*kleborate$resistance_score )
  coef_mat <- summary(model)$coefficients
  input<-c(colnames(fitness_mat)[i], coef_mat["kleborate$virulence_score", ],confint(model, "kleborate$virulence_score", level = 0.95) )
  names(input)<-c("Condition","Estimate","Std", "t_value","pval" , "Upper" , "Lower" )
  output<-rbind(output, input)
}
output<-data.frame(output)
output$Estimate <- as.numeric( as.character(output$Estimate) )
output$pval<- as.numeric( as.character(output$pval) )
output$Upper<- as.numeric( as.character(output$Upper) )
output$Lower<- as.numeric( as.character(output$Lower) )

output$p.adjusted <- p.adjust(output$pval, method = "BH")
row.names(output) <- NULL

output_short<-output[output$p.adjusted<0.05,] %>% arrange(Estimate)

ggplot(output_short, aes(x = reorder(Condition, Estimate), y = Estimate)) +
  geom_col(fill = "darkgreen", width = 0.7) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                width = 0.2, colour = "black", size = 0.8) +
  coord_flip() +
  labs(x = "Condition",
       y = "Fitness Effect") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8))


# Extract and plot the resistance-score coefficient from each per-condition model.
#resistance score 

output<-c()
for(i in 1:dim(fitness_mat)[2]){
  model<-lm( fitness_mat[,i] %>% pull(.) ~  kleborate$virulence_score + kleborate$resistance_score + kleborate$virulence_score*kleborate$resistance_score )
  coef_mat <- summary(model)$coefficients
  input<-c(colnames(fitness_mat)[i], coef_mat["kleborate$resistance_score", ],confint(model, "kleborate$resistance_score", level = 0.95) )
  names(input)<-c("Condition","Estimate","Std", "t_value","pval" , "Upper" , "Lower" )
  output<-rbind(output, input)
}
output<-data.frame(output)
output$Estimate <- as.numeric( as.character(output$Estimate) )
output$pval<- as.numeric( as.character(output$pval) )
output$Upper<- as.numeric( as.character(output$Upper) )
output$Lower<- as.numeric( as.character(output$Lower) )

output$p.adjusted <- p.adjust(output$pval, method = "BH")
row.names(output) <- NULL

output_short<-output[output$p.adjusted<0.05,] %>% arrange(Estimate)

ggplot(output_short, aes(x = reorder(Condition, Estimate), y = Estimate)) +
  geom_col(fill = "brown", width = 0.7) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                width = 0.2, colour = "black", size = 0.8) +
  coord_flip() +
  labs(x = "Condition",
       y = "Fitness Effect",
       title = "Resistance Score Correlation with Fitness Values") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# Extract and plot the virulence-by-resistance interaction coefficient.
#inetarction between resistance virulence 

output<-c()
for(i in 1:dim(fitness_mat)[2]){
  model<-lm( fitness_mat[,i] %>% pull(.) ~  kleborate$virulence_score + kleborate$resistance_score + kleborate$virulence_score*kleborate$resistance_score )
  coef_mat <- summary(model)$coefficients
  input<-c(colnames(fitness_mat)[i], coef_mat["kleborate$virulence_score:kleborate$resistance_score", ],confint(model, "kleborate$virulence_score:kleborate$resistance_score", level = 0.95) )
  names(input)<-c("Condition","Estimate","Std", "t_value","pval" , "Upper" , "Lower" )
  output<-rbind(output, input)
}
output<-data.frame(output)
output$Estimate <- as.numeric( as.character(output$Estimate) )
output$pval<- as.numeric( as.character(output$pval) )
output$Upper<- as.numeric( as.character(output$Upper) )
output$Lower<- as.numeric( as.character(output$Lower) )

output$p.adjusted <- p.adjust(output$pval, method = "BH")
row.names(output) <- NULL

output_short<-output[output$p.adjusted<0.05,] %>% arrange(Estimate)

ggplot(output_short, aes(x = reorder(Condition, Estimate), y = Estimate)) +
  geom_col(fill = "green", width = 0.7) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                width = 0.2, colour = "black", size = 0.8) +
  coord_flip() +
  labs(x = "Condition",
       y = "Fitness Effect") +
  theme_bw()  +
  theme(axis.text.y = element_text(size = 8))


# -----------------------------------------------------------------------------
# Sequence-type-level fitness summaries.
# Aggregate fitness by major STs, plot heatmaps, and correlate ST-level fitness
# with virulence/resistance summaries.
# -----------------------------------------------------------------------------
#STs
library(pheatmap)
fitness_mat<-read_csv("/Data/absolute_fitness_ML.csv")
# Approach 1: using sapply & var (numeric columns)
# Approach: check for >1 unique non-NA value (works generically)
colnames(fitness_mat)[1]<-"Name"


kleborate<-read_tsv("/Data/tot_results.tsv",show_col_types = FALSE)
kleborate$ST_short <-as.character(sapply(kleborate$ST, function(x) strsplit(x,"-")[[1]][1]))
kleborate$strain <- toupper(kleborate$strain)

shared<- intersect(fitness_mat$Name, kleborate$strain)
fitness_mat<-fitness_mat[match(shared, fitness_mat$Name), ]
kleborate<-kleborate[ match(shared, kleborate$strain),]


# Keep the most frequent STs and group all remaining STs as Other ST.
major_ST<-names(sort(-table(kleborate$ST_short)))[1:18]
kleborate$ST_short[which( !kleborate$ST_short %in%  major_ST)]<-"Other ST"


fitness_mat_combined <- fitness_mat
fitness_mat_combined$ST <- kleborate$ST_short

# 2. Summarise: group by ST, take mean of all numeric columns
# Average fitness values within each ST group.
fitness_summary <- fitness_mat_combined %>%
  group_by(ST) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

# 3. Prepare matrix for heatmap: remove the ST column and set rownames
mat <- fitness_summary %>%
  column_to_rownames(var = "ST") %>%
  as.matrix()

# 4. Optionally scale rows or columns if you want
# mat_scaled <- t(scale(t(mat)))  # scale rows

# 5. Plot heatmap

# Define breakpoints so 1 maps to the white colour
# Suppose your data ranges from minVal to maxVal
minVal <- min(mat, na.rm = TRUE)
maxVal <- max(mat, na.rm = TRUE)

# Create breaks that include 1 exactly
# For example, if your range is [0, 2]
breaks <- c(seq(minVal, 0.999, length.out = 50),
            1,
            seq(1.001, maxVal, length.out = 50))

# Create matching palette:
cols1   <- colorRampPalette(c("navy",  "white"))(50)   # below 1
cols2   <- colorRampPalette(c("white", "firebrick3"))(50)  # above 1
cols    <- c(cols1, cols2)



condition_tags<- read_csv("/Data/condition_group.csv", locale = locale(encoding = "UTF-8"),show_col_types = FALSE)
condition_tags$File_name <- gsub("[^A-Za-z0-9_]", "", condition_tags$File_name)
condition_tags$Clean_Tags <- iconv(condition_tags$Clean_Tags, from = "latin1", to = "UTF-8")

colnames(mat)  <- condition_tags$Clean_Tags[match(colnames(mat), condition_tags$File_name)]
mat<- mat[,!is.na(colnames(mat))]

Category <- condition_tags$Category[match(colnames(mat), condition_tags$Clean_Tags)]

#correlation




# Order columns by category
order_cols <- order(Category)

mat_ordered <- mat[, order_cols]
Category <- Category[order_cols]

summary_df <- kleborate %>%
  group_by(ST_short) %>%
  summarise(
    mean_score = mean(virulence_score, na.rm = TRUE),
    sd_score   = sd(virulence_score, na.rm = TRUE),
    n          = n(),
    se_score   = sd_score / sqrt(n),
    .groups = "drop"
  )


summary_df <- kleborate %>%
  group_by(ST_short) %>%
  summarise(
    mean_score = mean(resistance_score, na.rm = TRUE),
    sd_score   = sd(resistance_score, na.rm = TRUE),
    n          = n(),
    se_score   = sd_score / sqrt(n),
    .groups = "drop"
  )
#virulence_indicator<-c()
#for(i in 1:dim(mat_ordered)[2]){
#  virulence_indicator<-c(virulence_indicator, cor(summary_df$mean_score, mat_ordered[ match(summary_df$ST_short, row.names(mat_ordered) ) ,i] ))
#}


# Compute correlations + p-values
results <- data.frame(
  condition = colnames(mat_ordered),
  cor = NA,
  pval = NA
)

# Calculate Spearman correlation between ST-level score and each condition.
for (i in 1:ncol(mat_ordered)) {
  x <- summary_df$mean_score
  y <- mat_ordered[
    match(summary_df$ST_short, rownames(mat_ordered)),
    i
  ]
  
  test <- cor.test(x, y, method = "spearman")
  
  results$cor[i] <- test$estimate
  results$pval[i] <- test$p.value
}

# Keep original order
results$condition <- factor(results$condition, levels = colnames(mat_ordered))

# Significance
results <- results %>%
  mutate(significant = ifelse(pval < 0.05, "Significant", "Not significant"))

# Plot (vertical bars)
ggplot(results, aes(x = condition, y = cor, fill = significant)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Significant" = "red", "Not significant" = "steelblue")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Correlation with Virulence Score",
    x = "Condition",
    y = "Correlation"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7)
  )

results
# Convert table to data frame
results$Category <- Category
df <- as.data.frame(table(results$significant, results$Category))
colnames(df) <- c("Significance", "Category", "Count")

# Plot (100% stacked barplot)
ggplot(df, aes(x = Category, y = Count, fill = Significance)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("Significant" = "red", "Not significant" = "steelblue")) +
  labs(
    x = "Category",
    y = "Proportion"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

#correlation




# annotation for columns
annotation_col <- data.frame(Category = Category)
rownames(annotation_col) <- colnames(mat_ordered)

pheatmap(
  mat_ordered,
  cluster_rows = TRUE,
  cluster_cols = FALSE,   # disable x-axis clustering
  annotation_col = annotation_col,
  annotation_colors = category_colors,
  display_numbers = FALSE,
  color = cols,
  breaks = breaks,
  main = "Mean fitness by ST group",
  fontsize_row = 8,
  fontsize_col = 5
)

p <- pheatmap(
  mat_ordered,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = category_colors,
  display_numbers = FALSE,
  color = cols,
  breaks = breaks,
  main = "Mean fitness by ST group",
  fontsize_row = 8,
  fontsize_col = 8
)
row_order <- p$tree_row$order
ordered_rows <- rownames(mat_ordered)[row_order]
ordered_rows


kleborate$ST_short
kleborate$virulence_score
library(dplyr)
library(ggplot2)

summary_df <- kleborate %>%
  group_by(ST_short) %>%
  summarise(
    mean_score = mean(virulence_score, na.rm = TRUE),
    sd_score   = sd(virulence_score, na.rm = TRUE),
    n          = n(),
    se_score   = sd_score / sqrt(n),
    .groups = "drop"
  )

# enforce the heatmap row order
summary_df$ST_short <- factor(summary_df$ST_short, levels = ordered_rows )

ggplot(summary_df,
       aes(x = ST_short, y = mean_score)) +
  
  geom_col(fill = "yellow4", width = 0.7) +
  
  geom_errorbar(aes(ymin = mean_score - se_score,
                    ymax = mean_score + se_score),
                width = 0.2,
                linewidth = 0.6) +
  scale_x_discrete(limits = rev)+
  coord_flip() +
  
  labs(
    x = "Sequence Type (ST)",
    y = "Mean Virulence Score"
  ) +
  
  theme_bw(base_size = 14)



summary_df <- kleborate %>%
  group_by(ST_short) %>%
  summarise(
    mean_score = mean(resistance_score, na.rm = TRUE),
    sd_score   = sd(resistance_score, na.rm = TRUE),
    n          = n(),
    se_score   = sd_score / sqrt(n),
    .groups = "drop"
  )

# enforce the heatmap row order
summary_df$ST_short <- factor(summary_df$ST_short, levels = ordered_rows)

ggplot(summary_df,
       aes(x = ST_short, y = mean_score)) +
  
  geom_col(fill = "orange4", width = 0.7) +
  
  geom_errorbar(aes(ymin = mean_score - se_score,
                    ymax = mean_score + se_score),
                width = 0.2,
                linewidth = 0.6) +
  scale_x_discrete(limits = rev)+
  coord_flip() +
  
  labs(
    x = "Sequence Type (ST)",
    y = "Mean Resistance Score"
  ) +
  
  theme_bw(base_size = 14)


#Interactions
library(tidyverse)
# Directory containing IRIS colony-measurement output files.
dir_path <- "/Data/assay_results_full"
files <- list.files(dir_path,pattern = "-1_A\\.JPG\\.iris$")
files<-gsub("-1_A.JPG.iris","",files)
conditions<- unique(as.character( sapply(files, function(x) paste0( strsplit(x,"-")[[1]][1], "-",  strsplit(x,"-")[[1]][2]) ) ))  

cond_1_list<-c(37 ,37, 37, 37, 37,  37,  37,  37,  37,  61, 62,  70, 69)
cond_2_list<-c(49 ,50, 69, 70, 111, 112, 133, 191, 192, 166, 166, 62, 103)
cond_1_comb<-c(38 ,39, 40, 41, 42,  43,  44,  45,  46,  64, 65,  74, 75)


for(j in 1:length(cond_2_list)){
  print("------------")
  print(j)
  print(cond_1_list[j])
  print(conditions[cond_1_list[j]])
  
  print(cond_2_list[j])
  print(conditions[cond_2_list[j]])
  
  print(cond_1_comb[j])
  print(conditions[cond_1_comb[j]])
  print("------------")
}


#CeferdericolCeftazadime-0125,1ugml
#CeferdericolTigecycline-0625,2ugml
#ColistinMeropenem-02,006ugml
#ColistinMeropenem-04,006ugml


df_tot_comb<-c() 
for(j in 1:length(cond_2_list)){
  print(conditions[cond_1_comb[j]])
  LB<- read_csv(paste0(dir_path, "/", "condition_ICC_New_",conditions[127], ".csv"),show_col_types = FALSE)
  cond_1<-read_csv(paste0(dir_path, "/", "condition_ICC_New_",conditions[cond_1_list[j]], ".csv"))
  cond_2<-read_csv(paste0(dir_path, "/", "condition_ICC_New_",conditions[cond_2_list[j]], ".csv"))
  cond_comb<-read_csv(paste0(dir_path, "/", "condition_ICC_New_",conditions[cond_1_comb[j]], ".csv"))
  
  LB$label <- "Base"
  cond_1$label<-"Drug1"
  cond_2$label<-"Drug2"
  cond_comb$label<-"DrugComb"
  
  assay<-rbind(LB,cond_1,cond_2,cond_comb)
  
  # Rename columns properly
  LB2 <- LB %>%
    select(strain, median_fc_overall,label) %>%
    rename(fc_LB = median_fc_overall)
  
  cond_1_2 <- cond_1 %>%
    select(strain, median_fc_overall,label) %>%
    rename(fc_drug1 = median_fc_overall)
  
  cond_2_2 <- cond_2 %>%
    select(strain, median_fc_overall,label) %>%
    rename(fc_drug2 = median_fc_overall)
  
  cond_comb_2 <- cond_comb %>%
    select(strain, median_fc_overall,label) %>%
    rename(fc_comb = median_fc_overall)
  
  # Merge all four conditions
  df_all <- LB2 %>%
    inner_join(cond_1_2, by = "strain") %>%
    inner_join(cond_2_2, by = "strain") %>%
    inner_join(cond_comb_2, by = "strain")
  
  # Calculate expected phenotype under NO interaction
  eps <- 1e-6
  
  df_expected <- df_all %>%
    mutate(
      fc_expected = ((fc_drug1+eps) * (fc_drug2+eps) ) / (fc_LB+eps),
      interaction_score = fc_comb - fc_expected,
      interaction_log2 = log2( (fc_comb+eps) / fc_expected)
    )
  df_expected_short<- select(df_expected, strain, fc_expected)
  df_expected_short$label<-"DrugBoth"
  
  colnames(LB2) <- c("strain","fc","label")
  colnames(cond_1_2) <- c("strain","fc","label")
  colnames(cond_2_2) <- c("strain","fc","label")
  colnames(cond_comb_2) <- c("strain","fc","label")
  colnames(df_expected_short) <- c("strain","fc","label")
  
  df_tot<- rbind(LB2, cond_1_2,cond_2_2, cond_comb_2,df_expected_short)
  df_tot$condition<- conditions[cond_1_comb[j]]
  
  df_tot_comb<-rbind(df_tot_comb, df_tot ) 
}


condition_tags<- read_csv("/Data/condition_group.csv", locale = locale(encoding = "UTF-8"),show_col_types = FALSE)
condition_tags$File_name <- gsub("[^A-Za-z0-9_]", "", condition_tags$File_name)
condition_tags$Clean_Tags <- iconv(condition_tags$Clean_Tags, from = "latin1", to = "UTF-8")

df_tot_comb$condition <- gsub("[^A-Za-z0-9_]", "", df_tot_comb$condition)
df_tot_comb$Clean_Tags <- condition_tags$Clean_Tags[match(df_tot_comb$condition, condition_tags$File_name)]
df_tot_comb <-df_tot_comb[!is.na(df_tot_comb$Clean_Tags),]

ggplot(df_tot_comb,
       aes(x = Clean_Tags,
           y = fc,
           fill = label,
           linetype = label)) +
  
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7,
    outlier.size = 0.5
  ) +
  
  scale_fill_manual(values = c(
    "Base" = "grey95",
    "Drug1" = "#ADD8E6",
    "Drug2" = "#FFB6B6",
    "DrugBoth" = "#66C266",
    "DrugComb" = "white"
  )) +
  
  scale_linetype_manual(values = c(
    "Base" = "solid",
    "Drug1" = "solid",
    "Drug2" = "solid",
    "DrugBoth" = "solid",
    "DrugComb" = "solid"
  )) +
  
  theme_bw(base_size = 14) +
  ylim(0,3) +
  
  theme(
    axis.text.x = element_text(angle = 90, , vjust = 0.5, hjust = 1, size=12),
    panel.grid.minor = element_blank()
  ) +
  
  labs(
    x = "Condition",
    y = "Fitness (fc)",
    fill = "Group",
    linetype = "Group"
  )


df_tot_comb<-c() 
for(j in 1:length(cond_2_list)){
  print(conditions[cond_1_comb[j]])
  LB<- read_csv(paste0(dir_path, "/", "condition_ICC_New_",conditions[127], ".csv"),show_col_types = FALSE)
  cond_1<-read_csv(paste0(dir_path, "/", "condition_ICC_New_",conditions[cond_1_list[j]], ".csv"))
  cond_2<-read_csv(paste0(dir_path, "/", "condition_ICC_New_",conditions[cond_2_list[j]], ".csv"))
  cond_comb<-read_csv(paste0(dir_path, "/", "condition_ICC_New_",conditions[cond_1_comb[j]], ".csv"))
  
  LB$label <- "Base"
  cond_1$label<-"Drug1"
  cond_2$label<-"Drug2"
  cond_comb$label<-"DrugComb"
  
  # Rename columns properly
  LB2 <- LB %>%
    select(strain, median_fc_overall,label) %>%
    rename(fc_LB = median_fc_overall)
  
  cond_1_2 <- cond_1 %>%
    select(strain, median_fc_overall,label) %>%
    rename(fc_drug1 = median_fc_overall)
  
  cond_2_2 <- cond_2 %>%
    select(strain, median_fc_overall,label) %>%
    rename(fc_drug2 = median_fc_overall)
  
  cond_comb_2 <- cond_comb %>%
    select(strain, median_fc_overall,label) %>%
    rename(fc_comb = median_fc_overall)
  
  # Merge all four conditions
  df_all <- LB2 %>%
    inner_join(cond_1_2, by = "strain") %>%
    inner_join(cond_2_2, by = "strain") %>%
    inner_join(cond_comb_2, by = "strain")
  
  # Calculate expected phenotype under NO interaction
  eps <- 1e-6
  
  df_expected <- df_all %>%
    mutate(
      fc_expected = ((fc_drug1+eps) * (fc_drug2+eps) ) / (fc_LB+eps),
      interaction_score = fc_comb - fc_expected,
      interaction_log2 = log2( (fc_comb+eps) / fc_expected)
    )
  
  df_expected <- df_expected %>%
    mutate(
      log_comb = log2(fc_comb+eps),
      log_expected = log2(fc_expected+eps),
      interaction = log2((fc_comb+eps ) / (fc_expected +eps))
    )
  
  # robust estimate of noise using median absolute deviation
  sigma <- mad(df_expected$interaction, na.rm = TRUE)
  df_expected <- df_expected %>%
    mutate(
      z_score = interaction / sigma,
      p_value = 2 * pnorm(-abs(z_score))
    )
  
  df_expected <- df_expected %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH")
    )
  
  df_expected <- df_expected %>%
    mutate(
      interaction_type = case_when(
        
        p_adj < 0.05 & interaction < 0 ~ "Significant synergy",
        p_adj < 0.05 & interaction > 0 ~ "Significant antagonism",
        
        TRUE ~ "No significant interaction"
      )
    )
  
  print(wilcox.test(df_expected$log_comb, df_expected$log_expected, paired = TRUE)) 
  df_tmp<-data.frame(table(df_expected$interaction_type))
  df_tmp$condition<- conditions[cond_1_comb[j]]
  df_tot_comb<-rbind(df_tot_comb,df_tmp) 
}


colnames(df_tot_comb)<-c("Interaction_type", "Freq", "condition")
df_tot_comb


condition_tags<- read_csv("/Data/condition_group.csv", locale = locale(encoding = "UTF-8"),show_col_types = FALSE)
condition_tags$File_name <- gsub("[^A-Za-z0-9_]", "", condition_tags$File_name)
condition_tags$Clean_Tags <- iconv(condition_tags$Clean_Tags, from = "latin1", to = "UTF-8")

df_tot_comb$condition <- gsub("[^A-Za-z0-9_]", "", df_tot_comb$condition)
df_tot_comb$Clean_Tags <- condition_tags$Clean_Tags[match(df_tot_comb$condition, condition_tags$File_name)]
df_tot_comb <-df_tot_comb[!is.na(df_tot_comb$Clean_Tags),]



ggplot(df_tot_comb,
       aes(x = Clean_Tags,
           y = Freq,
           fill = Interaction_type)) +
  
  geom_col(width = 0.8, color = "black", linewidth = 0.2) +
  
  scale_fill_manual(values = c(
    "Significant synergy" = "#1B9E77",        # green
    "Significant antagonism" = "#D95F02",    # orange/red
    "No significant interaction" = "#BDBDBD" # grey
  )) +
  
  theme_bw(base_size = 14) +
  
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  
  labs(
    x = "Condition",
    y = "Number of interactions",
    fill = "Interaction type"
  )


#statisctal test 

results <- df_tot_comb %>%
  
  # keep only antagonism and synergy
  filter(Interaction_type %in% c("Significant antagonism",
                                 "Significant synergy")) %>%
  
  # convert to wide format
  pivot_wider(
    names_from = Interaction_type,
    values_from = Freq,
    values_fill = 0
  ) %>%
  
  # compute statistics per condition
  rowwise() %>%
  mutate(
    
    total = `Significant antagonism` + `Significant synergy`,
    
    enrichment = `Significant antagonism` / total,
    
    p_value = binom.test(
      x = `Significant antagonism`,
      n = total,
      p = 0.5,
      alternative = "greater"
    )$p.value
    
  ) %>%
  
  ungroup() %>%
  
  # adjust for multiple testing
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    
    dominant_interaction = case_when(
      enrichment > 0.5 & p_adj < 0.05 ~ "Antagonism enriched",
      enrichment < 0.5 & p_adj < 0.05 ~ "Synergy enriched",
      TRUE ~ "No significant difference"
    )
  ) %>%
  
  arrange(p_adj)

results$condition


#plasmid count 
plasmid<-read_tsv("/Data/summary_ARG.tsv")
plasmid$`#FILE` <- gsub("abricate_","", plasmid$`#FILE`)
plasmid$`#FILE` <- gsub(".fasta.tsv","", plasmid$`#FILE`)

Fitness<-read_csv("/Data/absolute_fitness_ML.csv")
colnames(Fitness)[1]<-"Name"
Fitness<-Fitness[match(plasmid$`#FILE`,Fitness$Name),]

#STs filtering
kleborate<-read_tsv("/Data/tot_results.tsv",show_col_types = FALSE)
kleborate$ST_short <-as.character(sapply(kleborate$ST, function(x) strsplit(x,"-")[[1]][1]))
ST_vector <- kleborate$ST_short[match(Fitness$Name, kleborate$strain)]
#Fitness<- Fitness[which(ST_vector=="ST2096"),]
#plasmid <- plasmid[which(ST_vector=="ST2096"),]

fitness_mat_combined <- Fitness
fitness_mat_combined$ST <- plasmid$NUM_FOUND

# 2. Summarise: group by ST, take mean of all numeric columns
fitness_summary <- fitness_mat_combined %>%
  group_by(ST) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

# 3. Prepare matrix for heatmap: remove the ST column and set rownames
mat <- fitness_summary %>%
  column_to_rownames(var = "ST") %>%
  as.matrix()

# 4. Optionally scale rows or columns if you want
# mat_scaled <- t(scale(t(mat)))  # scale rows

# 5. Plot heatmap

# Define breakpoints so 1 maps to the white colour
# Suppose your data ranges from minVal to maxVal
minVal <- min(mat, na.rm = TRUE)
maxVal <- max(mat, na.rm = TRUE)

# Create breaks that include 1 exactly
# For example, if your range is [0, 2]
breaks <- c(seq(minVal, 0.999, length.out = 50),
            1,
            seq(1.001, maxVal, length.out = 50))



condition_tags<- read_csv("/Data/condition_group.csv", locale = locale(encoding = "UTF-8"),show_col_types = FALSE)
condition_tags$File_name <- gsub("[^A-Za-z0-9_]", "", condition_tags$File_name)
condition_tags$Clean_Tags <- iconv(condition_tags$Clean_Tags, from = "latin1", to = "UTF-8")

colnames(correlation_matrix)


#mat <- as.matrix(tot_vals)
colnames(mat) <- condition_tags$Clean_Tags[match(colnames(mat), condition_tags$File_name)]
mat<-mat[,!is.na(colnames(mat))]
category <- condition_tags$Category[match(colnames(mat), condition_tags$Clean_Tags)]

# Create matching palette:
#cols1   <- colorRampPalette(c("navy",  "white"))(50)   # below 1
#cols2   <- colorRampPalette(c("white", "firebrick3"))(50)  # above 1
#cols    <- c(cols1, cols2)
# get ordering index
order_cat <- order(category)
# reorder matrix columns
mat_ordered <- mat[, order_cat]
# reorder category vector
category <- category[order_cat]

# annotation dataframe
annotation <- data.frame(Category = category)
rownames(annotation) <- colnames(mat_ordered)

# category colors
category_colors <- list(
  Category = c(
    "Chemical stress" = "#377EB8",
    "Carbon source utilisation" = "#4DAF4A",
    "Antibiotics" = "#E41A1C",
    "Environmental stress" = "#984EA3",
    "Antiseptic" = "#FF7F00",
    "Stress combination" = "#A65628",
    "Metal stress" = "#B8860B",
    "Base" = "#999999",
    "Envelope stress" = "#17BECF"
  )
)

# force white at zero
rg <- max(abs(mat_ordered))
breaks <- seq(-rg, rg, length.out = 101)

breaks <- c(seq(minVal, 0.999, length.out = 50),
            1,
            seq(1.001, maxVal, length.out = 50))

# Create matching palette:
cols1   <- colorRampPalette(c("navy",  "white"))(50)   # below 1
cols2   <- colorRampPalette(c("white", "firebrick3"))(50)  # above 1
cols    <- c(cols1, cols2)


# heatmap
pheatmap(
  mat_ordered,
  color          = cols,
  breaks         = breaks,
  annotation_col = annotation,
  annotation_colors = category_colors,
  cluster_rows = F,
  cluster_cols = F,
  scale = "none",
  fontsize_row = 8,
  fontsize_col = 6,
  border_color = NA,
  main = "Plasmid contribution to fitness"
)

#Summazry (Fix the issue)

# --------------------------------------------------
# 1. Normalize relative to 0 plasmids (row 1)
# --------------------------------------------------
mat_norm <- sweep(mat, 2, mat[1, ], "/")

# --------------------------------------------------
# 2. Convert to long format
# --------------------------------------------------
mat_long <- mat_norm %>%
  as.data.frame() %>%
  rownames_to_column("plasmid_count") %>%
  pivot_longer(
    cols = -plasmid_count,
    names_to = "condition",
    values_to = "value"
  ) %>%
  mutate(
    plasmid_count = as.numeric(plasmid_count)
  )

# --------------------------------------------------
# 3. Attach category (ASSUMES SAME ORDER AS COLS)
# --------------------------------------------------
#mat_long$category <- rep(category, each = nrow(mat_norm))
#condition_tags$Clean_Tags <- iconv(condition_tags$Clean_Tags, from = "latin1", to = "UTF-8")
mat_long$category<- condition_tags$Category[match(mat_long$condition, condition_tags$Clean_Tags)]

# --------------------------------------------------
# 4. Summarise per plasmid count + category
# --------------------------------------------------
summary_df <- mat_long %>%
  group_by(plasmid_count, category) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value   = sd(value, na.rm = TRUE),
    n          = n(),
    se         = sd_value / sqrt(n),
    ci_lower   = mean_value - 1.96 * se,
    ci_upper   = mean_value + 1.96 * se,
    .groups = "drop"
  )

# --------------------------------------------------
# 5. Colors
# --------------------------------------------------
category_colors <- c(
  "Chemical stress" = "#377EB8",
  "Carbon source utilisation" = "#4DAF4A",
  "Antibiotics" = "#E41A1C",
  "Environmental stress" = "#984EA3",
  "Antiseptic" = "#FF7F00",
  "Stress combination" = "#A65628",
  "Metal stress" = "#B8860B",
  "Base" = "#999999",
  "Envelope stress" = "#17BECF"
)

# --------------------------------------------------
# 6. Plot
# --------------------------------------------------
ggplot(summary_df,
       aes(x = plasmid_count,
           y = mean_value,
           color = category,
           fill = category)) +
  
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
              alpha = 0.2,
              color = NA) +
  
  geom_line(linewidth = 1.2) +
  
  geom_point(size = 2, alpha = 0.7, shape = 16) +
  
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +
  
  scale_color_manual(values = category_colors) +
  scale_fill_manual(values = category_colors) +
  
  scale_x_continuous(
    breaks = sort(unique(summary_df$plasmid_count))
  ) +
  
  labs(
    x = "Number of plasmids",
    y = "Relative fitness",
    color = "Category",
    fill = "Category"
  ) +
  
  theme_bw(base_size = 14)


ggplot(summary_df,
       aes(x = plasmid_count,
           y = mean_value,
           color = category,
           fill = category)) +
  
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
              alpha = 0.2,
              color = NA) +
  
  geom_line(linewidth = 1.2) +
  
  geom_point(size = 2, alpha = 0.7, shape = 16) +
  
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +
  ylim(range(0.5,3))+
  scale_color_manual(values = category_colors) +
  scale_fill_manual(values = category_colors) +
  
  scale_x_continuous(
    breaks = sort(unique(summary_df$plasmid_count))
  ) +
  
  labs(
    x = "Number of plasmids",
    y = "Relative fitness",
    color = "Category",
    fill = "Category"
  ) +
  
  theme_bw(base_size = 14)



#Feature count 






#Monday Schema 17

#17 - 18 Figures 

#18-25 write-up 







df_results
best_per_condition <- df_results %>%
  group_by(condition) %>%
  slice_max(order_by = rho_test, n = 1, with_ties = FALSE) %>%
  ungroup()

table(best_per_condition$model)
table(best_per_condition$mode)




boxplot(performance_df$rho_test ~ performance_df$model)
boxplot(performance_df$cov95_test ~ performance_df$model)


# Directory containing IRIS colony-measurement output files.
dir_path <- "/Data/assay_results_full"
files <- list.files(dir_path,pattern = "condition_ICC_New_")
ICC_vals<-c()
sd_vals<-c()

condition<-c()
for(i in 1:length(files)){
  tmp<-read_csv(paste0(dir_path, "/", files[i]),show_col_types = FALSE)
  tmp_tag<- gsub("condition_ICC_New_","",files[i])
  tmp_tag<- gsub(".csv","",tmp_tag)
  
  condition<-c(condition, gsub("[^A-Za-z0-9]", "", tmp_tag))
  ICC_vals<-c(ICC_vals, tmp$ICC_over_plates[1] )
  sd_vals<-c(sd_vals,sd(tmp$median_fc_overall))
}
ICC_df=data.frame(
  condition=condition,
  ICC_vals=ICC_vals,
  sd_vals=sd_vals
)

performance_df$ICC_val<- ICC_df$ICC_vals[match(performance_df$condition, ICC_df$condition)]
performance_df$rho_norm <- performance_df$rho_test / sqrt(performance_df$ICC_val)
hist(performance_df$rho_norm)

mean(ICC_vals[-1])

#write_csv(ICC_df, "/Data/ICC_df.csv" )


performance_df
best_models <- performance_df %>%
  group_by(condition) %>%
  slice_max(order_by = rho_norm, n = 1, with_ties = FALSE) %>%
  ungroup()

best_models

table(best_models$model)
ggplot(best_models, aes(x = reorder(condition, rho_test),
                        y = rho_test,
                        color = model,
                        shape = mode)) +
  geom_point(size = 3, alpha = 0.9) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Condition",
    y = expression(rho[test]),
    title = expression("Best model performance (" * rho[test] * ") per condition"),
    color = "Model",
    shape = "Mode"
  )


condition_group<- read_csv("/Data/condition_group.csv")
condition_group$File_name <- gsub("[^A-Za-z0-9_]", "", condition_group$File_name )

condition_group$Category[match(best_models$condition,condition_group$File_name)]
condition_group$File_name[match(best_models$condition,condition_group$File_name)]

best_models$Category<-condition_group$Category[match(best_models$condition,condition_group$File_name)]
table(best_models$Category)


library(ggplot2)

ggplot(best_models,
       aes(x = Category,
           y = rho_test,
           fill = Category)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(aes(color = Category),
              width = 0.2,
              alpha = 0.6,
              size = 2) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  labs(
    x = "Condition category",
    y = expression(rho[test]),
    title = expression("Prediction performance grouped by stress category")
  )


#Feature importance 


#write_csv(tot_importance, "/Data/SHAP_important_features.csv")
tot_importance <- read_csv("/Data/SHAP_important_features.csv")
tot_importance$feature_type <- ifelse(
  grepl("\\.\\.\\.\\.\\.GT", tot_importance$feature), "SNP",
  ifelse(
    grepl("^Clusters_", tot_importance$feature),
    "STR",
    "PAN"
  )
)


max_rank <- max(tot_importance$rank_final, na.rm = TRUE)


#conditions<-read_csv("/Data/condition_group.csv",show_col_types = FALSE)
#conditions$File_name <- gsub("[^A-Za-z0-9]+", "", conditions$File_name)
#tot_importance$category<-conditions$Category[match(tot_importance$target, conditions$File_name)]

condition_tags<- read_csv("/Data/condition_group.csv", locale = locale(encoding = "UTF-8"),show_col_types = FALSE)
condition_tags$File_name <- gsub("[^A-Za-z0-9_]", "", condition_tags$File_name)
condition_tags$Clean_Tags <- iconv(condition_tags$Clean_Tags, from = "latin1", to = "UTF-8")

tot_importance$category<- condition_tags$Category[match(tot_importance$target, condition_tags$File_name)]
tot_importance$Clean_Tags <- condition_tags$Clean_Tags[match(tot_importance$target, condition_tags$File_name)]

tot_importance <- tot_importance[!is.na(tot_importance$Clean_Tags),]

# Step 1: define category order
category_order <- unique(tot_importance$category)

# Step 2: create ordered target levels grouped by category
target_order <- tot_importance %>%
  distinct(target, category) %>%
  mutate(category = factor(category, levels = category_order)) %>%
  arrange(category, target) %>%   # group by category, then alphabetically
  pull(target)

# Step 3: apply ordering
tot_importance <- tot_importance %>%
  mutate(
    category = factor(category, levels = category_order),
    target   = factor(target, levels = target_order)
  )


category_colors <- list(
  Category = c(
    "Chemical stress" = "#377EB8",
    "Carbon source utilisation" = "#4DAF4A",
    "Antibiotics" = "#E41A1C",
    "Environmental stress" = "#984EA3",
    "Antiseptic" = "#FF7F00",
    "Stress combination" = "#A65628",
    "Metal stress" = "#B8860B",
    "Base" = "#999999",
    "Envelope stress" = "#17BECF"
  )
)

# attach category to each target
tot_importance <- tot_importance %>%
  left_join(condition_tags %>% select(Clean_Tags, Category),
            by = c("target" = "Clean_Tags"))

# order x-axis by category then target
tot_importance <- tot_importance %>%
  arrange(Category, target) %>%
  mutate(target = factor(target, levels = unique(target)))

# plot
ggplot(tot_importance,
       aes(x = target,
           y = rank_final,
           fill = feature_type)) +
  
  geom_tile(color = "white", linewidth = 0.3) +
  
  scale_y_reverse(
    breaks = 1:10,
    limits = c(max_rank + 0.5, 0.5),
    expand = c(0,0)
  ) +
  
  scale_fill_manual(
    values = c(
      "PAN" = "#66bd63",
      "SNP" = "#4575b4",
      "STR" = "#fdae61"
    ),
    drop = FALSE
  ) +
  
  facet_grid(. ~ category, scales = "free_x", space = "free_x") +
  
  labs(
    x = "Condition",
    y = "Feature rank",
    fill = "Feature type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
    panel.spacing.x = unit(0, "pt"),
    panel.spacing.y = unit(0, "pt"),
    strip.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)
  )


library(dplyr)
library(ggplot2)

# reorder x-axis by category
tot_importance <- tot_importance %>%
  arrange(category, target) %>%
  mutate(target = factor(target, levels = unique(target)))

# category band data
cat_band <- tot_importance %>%
  distinct(target, category) %>%
  mutate(rank_final = 0)

# plot
ggplot(tot_importance, aes(x = target, y = rank_final)) +
  
  # heatmap tiles
  geom_tile(aes(fill = feature_type),
            color = "white",
            linewidth = 0.3) +
  
  # category band on top
  geom_tile(data = cat_band,
            aes(x = target, y = 0, fill = category),
            height = 0.8) +
  
  scale_y_reverse(
    breaks = 1:10,
    limits = c(max_rank + 0.5, -0.5),
    expand = c(0,0)
  ) +
  
  scale_fill_manual(
    values = c(
      category_colors$Category,
      "PAN" = "#66bd63",
      "SNP" = "#4575b4",
      "STR" = "#fdae61"
    )
  ) +
  
  labs(
    x = "Condition",
    y = "Feature rank",
    fill = ""
  ) +
  theme_bw(base_size = 14)+
  #theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
    panel.grid = element_blank()
  )

levels(tot_importance$target)
colnames(tot_importance)
#write_csv(tot_importance, "/Data/predictive_totall_importance_1.csv")

x_order <- tot_importance %>%
  arrange(category, target) %>%
  distinct(target) %>%
  pull(target)

#features 
inp<-read_csv("/Data/predictive_features.csv")
colnames(inp)[1]<-"Name"

condition_tags<- read_csv("/Data/condition_group.csv", locale = locale(encoding = "UTF-8"),show_col_types = FALSE)
condition_tags$File_name <- gsub("[^A-Za-z0-9_]", "", condition_tags$File_name)
condition_tags$Clean_Tags <- iconv(condition_tags$Clean_Tags, from = "latin1", to = "UTF-8")


inp$category<- condition_tags$Category[match(inp$Name, condition_tags$File_name)]
inp$Clean_Tags <- condition_tags$Clean_Tags[match(inp$Name, condition_tags$File_name)]
inp <- inp[!is.na(inp$Clean_Tags),]
#write_csv(tot_importance, "/Data/predictive_totall_importance.csv")



# reuse the exact order from the heatmap
inp_plot <- inp %>%
  mutate(Name = factor(Name, levels = x_order))

# plot
ggplot(inp_plot, aes(x = Name, y = Pangenome, fill = category)) +
  
  geom_bar(stat = "identity", width = 0.8) +
  
  scale_fill_manual(values = category_colors$Category) +
  
  labs(
    x = "Condition",
    y = "Pangenome size",
    fill = "Category"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5), 
    panel.grid.minor = element_blank()
  )

ggplot(inp_plot, aes(x = Name, y = SNP, fill = category)) +
  
  geom_bar(stat = "identity", width = 0.8) +
  
  scale_fill_manual(values = category_colors$Category) +
  
  labs(
    x = "Condition",
    y = "SNP size",
    fill = "Category"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5), 
    panel.grid.minor = element_blank()
  )

#Virulence and resistance score 
library(tidyverse)
ARG <- read_csv("/Data/AMR_genes.csv")
colnames(ARG)[1]<-"Name"

kleborate<-read_tsv("/Data/tot_results.tsv",show_col_types = FALSE)
kleborate<-kleborate[match( ARG$Name, kleborate$strain),]

Fitness<-read_csv("/Data/absolute_fitness_mat_ML_22Feb.csv")

colnames(Fitness)[1]<-"Name"
Fitness<-Fitness[match(ARG$Name,Fitness$Name),]

kleborate$combined <- paste0("VirScore=",kleborate$virulence_score, "_", "ResScore=",kleborate$resistance_score)
Fitness$combined <- kleborate$combined
Fitness<-Fitness[Fitness$combined != "VirScore=NA_ResScore=NA" & Fitness$combined != "VirScore=5_ResScore=2" & Fitness$combined != "VirScore=2_ResScore=0",]

condition_tags<- read_csv("/Data/condition_group.csv", locale = locale(encoding = "UTF-8"),show_col_types = FALSE)
condition_tags$File_name <- gsub("[^A-Za-z0-9_]", "", condition_tags$File_name)
condition_tags$Clean_Tags <- iconv(condition_tags$Clean_Tags, from = "latin1", to = "UTF-8")

colnames(Fitness) <- condition_tags$Clean_Tags[match(colnames(Fitness), condition_tags$File_name)]
colnames(Fitness)[length(colnames(Fitness))] <- "combined"
Fitness<- Fitness[,!is.na(colnames(Fitness))]


inp$Clean_Tags <- condition_tags$Clean_Tags[match(inp$Name, condition_tags$File_name)]

#inp$Clean_Tags <- condition_tags$Clean_Tags[match(inp$Name, condition_tags$File_name)]
#inp <- inp[!is.na(inp$Clean_Tags),]

#df<-as.data.frame.matrix(table(Fitness$combined, Fitness$ST_short)) 
#write.csv(df, "/Data/df_ST.csv")

fitness_mat_combined <- Fitness
fitness_mat_combined$ST <- Fitness$combined

# 2. Summarise: group by ST, take mean of all numeric columns
fitness_summary <- fitness_mat_combined %>%
  group_by(ST) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

# 3. Prepare matrix for heatmap: remove the ST column and set rownames
mat <- fitness_summary %>%
  column_to_rownames(var = "ST") %>%
  as.matrix()


Category<- condition_tags$Category[match(colnames(mat), condition_tags$Clean_Tags)]

annotation_col <- data.frame(
  Category = Category
)
rownames(annotation_col) <- colnames(mat)

# define desired category order (optional but recommended)
category_order <- c(
  "Chemical stress",
  "Carbon source utilisation",
  "Antibiotics",
  "Environmental stress",
  "Antiseptic",
  "Stress combination",
  "Metal stress",
  "Base",
  "Envelope stress"
)

annotation_col$Category <- factor(annotation_col$Category,
                                  levels = category_order)

# reorder matrix columns
mat <- mat[, order(annotation_col$Category)]
annotation_col <- annotation_col[colnames(mat), , drop = FALSE]

annotation_colors <- list(
  Category = category_colors
)

# 4. Optionally scale rows or columns if you want
# mat_scaled <- t(scale(t(mat)))  # scale rows

# 5. Plot heatmap

# Define breakpoints so 1 maps to the white colour
# Suppose your data ranges from minVal to maxVal
minVal <- min(mat, na.rm = TRUE)
maxVal <- max(mat, na.rm = TRUE)

breaks <- c(seq(minVal, 0.999, length.out = 50),
            1,
            seq(1.001, maxVal, length.out = 50))

# fix above code 
Category <- condition_tags$Category[
  match(colnames(mat), condition_tags$Clean_Tags)
]

annotation_col <- data.frame(Category = Category)
rownames(annotation_col) <- colnames(mat)

# --------------------------------------------------
# 2. Clean categories
# --------------------------------------------------
annotation_col$Category <- trimws(annotation_col$Category)

# remove columns with missing category
valid <- !is.na(annotation_col$Category)
mat <- mat[, valid]
annotation_col <- annotation_col[valid, , drop = FALSE]

# --------------------------------------------------
# 3. Set category order
# --------------------------------------------------
category_order <- c(
  "Chemical stress",
  "Carbon source utilisation",
  "Antibiotics",
  "Environmental stress",
  "Antiseptic",
  "Stress combination",
  "Metal stress",
  "Base",
  "Envelope stress"
)

annotation_col$Category <- factor(
  annotation_col$Category,
  levels = category_order
)

# reorder columns
mat <- mat[, order(annotation_col$Category)]
annotation_col <- annotation_col[colnames(mat), , drop = FALSE]

# --------------------------------------------------
# 4. Define colors (correct structure)
# --------------------------------------------------
category_colors <- list(
  Category = c(
    "Chemical stress" = "#377EB8",
    "Carbon source utilisation" = "#4DAF4A",
    "Antibiotics" = "#E41A1C",
    "Environmental stress" = "#984EA3",
    "Antiseptic" = "#FF7F00",
    "Stress combination" = "#A65628",
    "Metal stress" = "#B8860B",
    "Base" = "#999999",
    "Envelope stress" = "#17BECF"
  )
)

# --------------------------------------------------
# 5. Color scale centered at 1
# --------------------------------------------------
minVal <- min(mat, na.rm = TRUE)
maxVal <- max(mat, na.rm = TRUE)

breaks <- c(
  seq(minVal, 0.999, length.out = 50),
  1,
  seq(1.001, maxVal, length.out = 50)
)

cols <- colorRampPalette(c("blue", "white", "red"))(100)

# --------------------------------------------------
# 6. Plot heatmap
# --------------------------------------------------
pheatmap(
  mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = category_colors,
  color = cols,
  breaks = breaks,
  main = "Mean fitness by condition",
  fontsize_row = 8,
  fontsize_col = 6,
  border_color = NA
)




# For example, if your range is [0, 2]


library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(tibble)

category_colors <- c(
  "Chemical stress" = "#377EB8",
  "Carbon source utilisation" = "#4DAF4A",
  "Antibiotics" = "#E41A1C",
  "Environmental stress" = "#984EA3",
  "Antiseptic" = "#FF7F00",
  "Stress combination" = "#A65628",
  "Metal stress" = "#B8860B",
  "Base" = "#999999",
  "Envelope stress" = "#17BECF"
)

mat %>%
  as.data.frame() %>%
  rownames_to_column("ST") %>%
  pivot_longer(-ST, names_to = "Clean_Tags", values_to = "fitness") %>%
  mutate(
    Category = condition_tags$Category[
      match(Clean_Tags, condition_tags$Clean_Tags)
    ],
    Virulence = as.numeric(str_extract(ST, "(?<=VirScore=)\\d+")),
    Resistance = as.numeric(str_extract(ST, "(?<=ResScore=)\\d+"))
  ) %>%
  group_by(Category, Virulence, Resistance) %>%
  summarise(mean_fitness = mean(fitness, na.rm = TRUE),
            .groups = "drop") %>%
  ggplot(aes(x = Virulence,
             y = Resistance,
             size = mean_fitness,
             fill = Category)) +
  
  geom_point(shape = 21, color = "black", alpha = 0.85) +
  
  scale_fill_manual(values = category_colors) +
  
  scale_size_continuous(range = c(2, 10)) +
  
  facet_wrap(~Category) +
  
  theme_bw(base_size = 14) +
  
  theme(
    panel.grid.minor = element_blank()
  ) +
  
  labs(
    x = "Virulence score",
    y = "Resistance score",
    size = "Mean fitness",
    title = "Fitness landscape across virulence and resistance scores by category"
  )

#Pleiotropy SNP and Gene 

#pleiotropy 
gwas_hits <- read_csv("/Data/GWAS_results_SNPs_Pangenome.csv")

condition_tags<- read_csv("/Data/condition_group.csv", locale = locale(encoding = "UTF-8"),show_col_types = FALSE)
condition_tags$File_name <- gsub("[^A-Za-z0-9_]", "", condition_tags$File_name)
condition_tags$Clean_Tags <- iconv(condition_tags$Clean_Tags, from = "latin1", to = "UTF-8")

#colnames(Fitness) <- condition_tags$Clean_Tags[match(colnames(Fitness), condition_tags$File_name)]
#colnames(Fitness)[length(colnames(Fitness))] <- "combined"
#Fitness<- Fitness[,!is.na(colnames(Fitness))]

#conditions<-read_csv("/Data/condition_group.csv",show_col_types = FALSE)
#conditions$File_name <- gsub("[^A-Za-z0-9]+", "", conditions$File_name)

gwas_hits$category<-condition_tags$Category[match(gwas_hits$condition, condition_tags$File_name)]
gwas_hits$condition<-condition_tags$Clean_Tags[match(gwas_hits$condition, condition_tags$File_name)]

gwas_hits <- gwas_hits[ !is.na(gwas_hits$condition) ,]

gwas_hits<- gwas_hits[gwas_hits$type=="Pangenome",]
gwas_hits<- select(gwas_hits,Gene,Odds_ratio,`Non-unique Gene name`,Annotation ,condition, category )



pleiotropy_df <- gwas_hits %>%
  group_by(Gene, `Non-unique Gene name`,Annotation) %>%
  summarise(
    n_conditions = n_distinct(condition),
    n_categories = n_distinct(category),
    conditions = paste(unique(condition), collapse = ";"),
    categories = paste(unique(category), collapse = ";"),
    pleiotropy_score = n_categories,  # main score
    .groups = "drop"
  )

pleiotropic_genes <- pleiotropy_df %>%
  filter(n_conditions > 1)

pleiotropic_category_genes <- pleiotropy_df %>%
  filter(n_categories > 1)

max_categories <- n_distinct(gwas_hits$category)

pleiotropy_df <- pleiotropy_df %>%
  mutate(
    pleiotropy_score_norm = n_categories / max_categories
  )

#top_genes <- pleiotropy_df %>%
#  arrange(desc(pleiotropy_score), desc(n_conditions)) %>%
#  head(20)

top_genes <- pleiotropy_df %>%
  arrange(desc(pleiotropy_score), desc(n_conditions)) 


top_genes <- top_genes[top_genes$pleiotropy_score==9,]
write_csv(top_genes, "/Data/top_pleiotropy_genes.csv")


gwas_hits_short<- gwas_hits[which(gwas_hits$Gene %in% top_genes$Gene),]

gwas_hits_short<- select(gwas_hits_short, Gene,Odds_ratio, condition,category )

library(dplyr)
library(ggplot2)

# -------------------------
# Prepare data (ORDER FIRST)
# -------------------------

plot_df <- gene_pleiotropy_avg %>%
  arrange(category, label) %>%
  mutate(
    label = factor(label, levels = unique(label)),
    x = as.numeric(label)   # 🔴 numeric positions for vline
  )

# -------------------------
# Compute boundaries
# -------------------------

boundaries <- plot_df %>%
  distinct(label, category, x) %>%
  group_by(category) %>%
  summarise(max_x = max(x), .groups = "drop") %>%
  pull(max_x)

# -------------------------
# Plot
# -------------------------

ggplot(plot_df,
       aes(x = label,
           y = LOCUS_TAG,
           fill = mean_Odds_ratio)) +
  
  geom_tile(color = "white", linewidth = 0.1) +
  
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 1,
    na.value = "grey85",
    name = "Odds ratio"
  ) +
  
  # 🔴 THIS NOW WORKS
  geom_vline(
    xintercept = boundaries[-length(boundaries)] + 0.5,
    color = "black",
    linewidth = 0.4
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank()
  ) +
  
  labs(
    x = "Condition",
    y = "Locus tag",
    title = "Gene pleiotropy heatmap"
  )
library(dplyr)
library(ggplot2)

# =========================
# DATA PREP
# =========================

category_order <- c(
  "Chemical stress",
  "Carbon source utilisation",
  "Antibiotics",
  "Environmental stress",
  "Antiseptic",
  "Stress combination",
  "Metal stress",
  "Base",
  "Envelope stress"
)

gwas_df <- gwas_hits_short %>%
  mutate(
    category = factor(category, levels = category_order)
  ) %>%
  arrange(category, condition) %>%
  mutate(
    condition = factor(condition, levels = unique(condition)),
    x = as.numeric(condition)   # 🔴 critical
  )

# =========================
# CATEGORY BOUNDARIES
# =========================

boundaries <- gwas_df %>%
  distinct(condition, category, x) %>%
  group_by(category) %>%
  summarise(max_x = max(x), .groups = "drop") %>%
  pull(max_x)

# =========================
# PLOT
# =========================

ggplot(gwas_df,
       aes(x = condition,
           y = Gene,
           fill = Odds_ratio)) +
  
  geom_tile(color = "white", linewidth = 0.1) +
  
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "grey",
    high = "#B2182B",
    midpoint = 1,
    na.value = "white",   # ✅ NA FIX
    name = "Odds ratio"
  ) +
  
  # ✅ separators now work
  geom_vline(
    xintercept = boundaries[-length(boundaries)] + 0.5,
    color = "black",
    linewidth = 0.4
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=5.8),
    panel.grid = element_blank()
  ) +
  
  labs(
    x = "Condition",
    y = "Gene"
  )

pangenome <- read_csv("/Data/gene_presence_absence_roary.csv")
gwas_df$Gene_Annotation <- pangenome$Annotation[match( gwas_df$Gene, pangenome$Gene)]

df<-cbind(unique(gwas_df$Gene), pangenome$Annotation[match( unique(gwas_df$Gene), pangenome$Gene)]) 
data.frame(df)

ggplot(gwas_df,
       aes(x = condition,
           y = Gene_Annotation,
           fill = Odds_ratio)) +
  
  geom_tile(color = "white", linewidth = 0.1) +
  
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "grey",
    high = "#B2182B",
    midpoint = 1,
    na.value = "white",   # ✅ NA FIX
    name = "Odds ratio"
  ) +
  
  # ✅ separators now work
  geom_vline(
    xintercept = boundaries[-length(boundaries)] + 0.5,
    color = "black",
    linewidth = 0.4
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=5.8),
    panel.grid = element_blank()
  ) +
  
  labs(
    x = "Condition",
    y = "Gene"
  )

unique(gwas_df$Gene)
#SNP Pleiotropy 
library(tidyverse)
library(clusterProfiler)
library(KEGGREST)

input <- read_csv("/Data/SNP_GWAS_annot.csv")
input$label <- gsub("[^A-Za-z0-9]+", "", input$label)

condition_tags<- read_csv("/Data/condition_group.csv", locale = locale(encoding = "UTF-8"),show_col_types = FALSE)
condition_tags$File_name <- gsub("[^A-Za-z0-9_]", "", condition_tags$File_name)
condition_tags$Clean_Tags <- iconv(condition_tags$Clean_Tags, from = "latin1", to = "UTF-8")

#condition_group<- read_csv("/Data/condition_group.csv")

#condition_group$File_name <- gsub("[^A-Za-z0-9]+", "", condition_group$File_name)

input$category<- condition_tags$Category[match(input$label, condition_tags$File_name)]
input <- input[!is.na(input$category),]

input_short<- select(input,label, Gene, LOCUS_TAG,EFFECT,PRODUCT, category,Odds_ratio)
head(input_short)

input_short$PRODUCT


gene_pleiotropy <- input_short %>%
  group_by(LOCUS_TAG,PRODUCT) %>%
  summarise(
    n_conditions = n_distinct(label),
    n_categories = n_distinct(category),
    Odds_ratio = mean(Odds_ratio),
    conditions = paste(unique(label), collapse = ";"),
    categories = paste(unique(category), collapse = ";"),
    pleiotropy_score = n_categories,
    .groups = "drop"
  ) %>%
  arrange(desc(pleiotropy_score))


arrange(gene_pleiotropy, pleiotropy_score)

top_genes <- gene_pleiotropy[gene_pleiotropy$pleiotropy_score==9,]
#write_csv(top_genes, "/Data/SNP_pleiotroy.csv")


top_genes <- gene_pleiotropy %>%
  arrange(desc(pleiotropy_score)) %>%
  head(100) %>%
  pull(LOCUS_TAG)

top_genes <- gene_pleiotropy$LOCUS_TAG[gene_pleiotropy$pleiotropy_score==9]


gene_pleiotropy_odds <- input_short[which( input_short$LOCUS_TAG %in% top_genes),]

gene_pleiotropy_avg <- gene_pleiotropy_odds %>%
  mutate(
    Odds_ratio = ifelse(is.infinite(Odds_ratio), NA, Odds_ratio)
  ) %>%
  group_by(LOCUS_TAG,PRODUCT,label,category) %>%
  summarise(
    mean_Odds_ratio = mean(Odds_ratio, na.rm = TRUE),
    n = n(),
    n_non_missing = sum(!is.na(Odds_ratio)),
    .groups = "drop"
  )

library(dplyr)
library(tidyr)
library(ggplot2)

category_order <- c(
  "Chemical stress",
  "Carbon source utilisation",
  "Antibiotics",
  "Environmental stress",
  "Antiseptic",
  "Stress combination",
  "Metal stress",
  "Base",
  "Envelope stress"
)

plot_df <- gene_pleiotropy_avg %>%
  arrange(category, label) %>%
  mutate(
    label = factor(label, levels = unique(label)),
    x = as.numeric(label)   # 🔴 numeric positions for vline
  )

# -------------------------
# Compute boundaries
# -------------------------

boundaries <- plot_df %>%
  distinct(label, category, x) %>%
  group_by(category) %>%
  summarise(max_x = max(x), .groups = "drop") %>%
  pull(max_x)

# -------------------------
# Plot
# -------------------------

ggplot(plot_df,
       aes(x = label,
           y = LOCUS_TAG,
           fill = mean_Odds_ratio)) +
  
  geom_tile(color = "white", linewidth = 0.1) +
  
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "grey",
    high = "#B2182B",
    midpoint = 1,
    na.value = "white",
    name = "Odds ratio"
  ) +
  
  # 🔴 THIS NOW WORKS
  geom_vline(
    xintercept = boundaries[-length(boundaries)] + 0.5,
    color = "black",
    linewidth = 0.4
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    axis.text.y = element_text( size = 4.5),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5),
    panel.grid = element_blank()
  ) +
  
  labs(
    x = "Condition",
    y = "Locus tag",
    title = "Gene pleiotropy heatmap"
  )


kegg_map <- read_tsv("/Data/kegg_map.tsv", col_names = c("LOCUS_TAG", "KEGG_ID", "Pathway"))
top_kegg <- kegg_map %>%
  filter(LOCUS_TAG %in% top_genes) %>%  # only your input genes
  pull(KEGG_ID) %>%
  unique()


# All background gene→pathway pairs
bg <- kegg_map %>%
  distinct(LOCUS_TAG, KEGG_ID)

# Count of all genes
N <- n_distinct(bg$LOCUS_TAG)

# Count of top genes
top_set <- kegg_map %>%
  filter(LOCUS_TAG %in% top_genes) %>%
  distinct(LOCUS_TAG) %>%
  pull()

n_top <- length(top_set)

# Prepare pathway counts
pathway_counts <- bg %>%
  group_by(KEGG_ID) %>%
  summarise(total_genes = n_distinct(LOCUS_TAG))

# In top set
top_pathway_counts <- kegg_map %>%
  filter(LOCUS_TAG %in% top_set) %>%
  group_by(KEGG_ID) %>%
  summarise(top_genes = n_distinct(LOCUS_TAG))

# Merge & fill
enrich_df <- pathway_counts %>%
  left_join(top_pathway_counts, by = "KEGG_ID") %>%
  mutate(top_genes = replace_na(top_genes, 0))

# Hypergeometric test
enrich_df <- enrich_df %>%
  rowwise() %>%
  mutate(
    p_value = phyper(top_genes - 1, total_genes, N - total_genes, n_top, lower.tail = FALSE),
    pathway = kegg_map$Pathway[match(KEGG_ID, kegg_map$KEGG_ID)]
  )

# Adjust for multiple testing (FDR)
enrich_df <- enrich_df %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

# View results
enrich_df


#GO enrichment 
go_map <- read_tsv("/Data/go_map.tsv", col_names = c("LOCUS_TAG", "GO", "Description"), comment = "#")

bg_go <- go_map %>%
  distinct(LOCUS_TAG, GO)

top_go <- go_map %>%
  filter(LOCUS_TAG %in% top_genes) %>% 
  distinct(LOCUS_TAG, GO)
# Count total background genes
N <- n_distinct(bg_go$LOCUS_TAG)

# Count total top (pleiotropic) genes
n_top <- n_distinct(top_go$LOCUS_TAG)

# Gene→GO frequency in background
bg_counts <- bg_go %>%
  group_by(GO) %>%
  summarise(bg_gene_count = n())

# Gene→GO frequency in target
top_counts <- top_go %>%
  group_by(GO) %>%
  summarise(top_gene_count = n())

# Combine and run hypergeometric test
enrich_go <- bg_counts %>%
  left_join(top_counts, by = "GO") %>%
  replace_na(list(top_gene_count = 0)) %>%
  rowwise() %>%
  mutate(
    p_value = phyper(
      top_gene_count - 1,
      bg_gene_count,
      N - bg_gene_count,
      n_top,
      lower.tail = FALSE
    )
  ) %>%
  ungroup()

# Adjust for multiple testing (FDR)
enrich_go <- enrich_go %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

enrich_go



snp_pleiotropy <- input_short %>%
  group_by(label,PRODUCT,LOCUS_TAG) %>%
  summarise(
    n_genes = n_distinct(Gene),
    n_conditions = n_distinct(label),      # same as 1 for each label
    n_categories = n_distinct(category),
    conditions = paste(unique(label), collapse = ";"),
    categories = paste(unique(category), collapse = ";"),
    pleiotropy_score = n_categories,
    .groups = "drop"
  )%>%
  arrange(desc(pleiotropy_score))

snp_pleiotropy


#plasmid character
kleborate<-read_tsv("/Data/tot_results.tsv",show_col_types = FALSE)

plasmid<-read_tsv("/Data/summary_ARG.tsv")
plasmid$`#FILE` <- gsub("abricate_","", plasmid$`#FILE`)
plasmid$`#FILE` <- gsub(".fasta.tsv","", plasmid$`#FILE`)

Fitness<-read_csv("/Data/absolute_fitness_ML.csv")
colnames(Fitness)[1]<-"Name"
Fitness<-Fitness[match(plasmid$`#FILE`,Fitness$Name),]

kleborate$strain<- toupper(kleborate$strain)
kleborate_short <-kleborate[ match(Fitness$Name,kleborate$strain),]

plasmid_short <- plasmid[match(Fitness$Name,plasmid$`#FILE`),]
colnames(plasmid_short)[1] <- "Name"

kleborate_short$strain

plasmid_rep <- c("Col440I_1","ColKP3_1", "ColRNAI_1", "IncFIA(HI1)_1_HI1","IncFIB(K)_1_Kpn3", "IncFIB(Mar)_1_pNDM-Mar", "IncFII_1_pKP91","IncHI1B_1_pNDM-MAR", "IncL/M(pOXA-48)_1_pOXA-48", "IncR_1"  )
carb_genes <- c("NDM-5","OXA-232","KPC-2", "OXA-48", "NDM-1","OXA-181")

for(i in 1:length(plasmid_rep)){
  for(j in 1:length(carb_genes)){
    tmp <- kleborate_short$Bla_Carb_acquired
    tmp[which(grepl(carb_genes[j], tmp))] <- 1
    tmp[which(tmp != 1 )] <- 0
    tmp<- as.numeric(tmp)
    
    tmp1<-plasmid_short[, which(colnames(plasmid_short)==plasmid_rep[i])] %>% pull(.)
    tmp1 <- ifelse(tmp1==".",0,1)
    
    #tmp1<-plasmid_short[, which(colnames(plasmid_short)=="IncL/M(pOXA-48)_1_pOXA-48")] %>% pull(.)
    
    
    print("-----------------")
    print(carb_genes[j])
    print(plasmid_rep[i])
    print(cor.test(tmp, tmp1))
    print("-----------------")
  }
}



#28-29 Introduction-Methods 
#30-4 Results 
#5-6 Discussion 

#3.45 4.46 GBP
#3.56 

#investment products
#portfolio 
#passive 
#strcutured notes 

#Global  #travel 
#meet in london
#permitted 

#strcutured note 150K USD 250K GBP

#direct exposure 
#Global investment into indices 


#2 April Introduction 
#3-10 Results/Discussion 
#10-15 Review 

#16-20 EHR 


#SNP GWAS
#filtering sNPs in R
dir_path <- "/Data/SNP_GWAS"

files <- list.files(
  path = dir_path,
  pattern = ".results.csv" , 
  recursive = TRUE,
  full.names = TRUE
)

files_short <- list.files(
  path = dir_path,
  pattern = ".results.csv" , 
  recursive = TRUE
)

files_short <- as.character( sapply( files_short, function(x) strsplit(x, "\\/")[[1]][2] ) ) 
files_short <- gsub(".results.csv","",files_short)

significant<-c()
for(i in 1:length(files_short)){
  print(i/length(files_short))
  tmp<- read_csv(files[i],show_col_types = FALSE)
  tmp$condition<-files_short[i]
 # tmp_1<-tmp[tmp$Benjamini_H_p < 0.05 & tmp$Best_pairwise_comp_p >= 0.05 & tmp$Worst_pairwise_comp_p >= 0.05,]
  tmp_2<-tmp[tmp$Benjamini_H_p < 0.05 & tmp$Best_pairwise_comp_p < 0.05 & tmp$Worst_pairwise_comp_p >= 0.05,]
 # tmp_3<-tmp[tmp$Benjamini_H_p < 0.05 & tmp$Best_pairwise_comp_p >= 0.05 & tmp$Worst_pairwise_comp_p < 0.05,]
  tmp_2_unique <- unique(paste(tmp_2$Number_pos_present_in,tmp_2$Number_neg_present_in) )
  
  tmp_4<-tmp[tmp$Benjamini_H_p < 0.05 & tmp$Best_pairwise_comp_p < 0.05 & tmp$Worst_pairwise_comp_p < 0.05,]
  tmp_4_unique <- unique(paste(tmp_4$Number_pos_present_in,tmp_4$Number_neg_present_in) )
  #significant<-rbind(significant,c( "0", dim(tmp)[1],files_short[i] ))
  #significant<-rbind(significant,c( "1", dim(tmp_1)[1],files_short[i] ))
  significant<-rbind(significant,c( "Broad",length(tmp_2_unique),files_short[i] ))
  #significant<-rbind(significant,c( "3", dim(tmp_3)[1],files_short[i] ))
  significant<-rbind(significant,c( "Narrow",length(tmp_4_unique),files_short[i] ))
}
significant<- data.frame(significant)
colnames(significant) <- c("Type","Count","Condition")

significant_filtered<- significant[significant$Type == "Broad",]

significant_filtered$Condition <- gsub("[^A-Za-z0-9_]", "", significant_filtered$Condition )
significant_filtered$Condition <- gsub("_", "", significant_filtered$Condition )

condition_tags<- read_csv("/Data/condition_group.csv", locale = locale(encoding = "UTF-8"))
condition_tags$File_name <- gsub("[^A-Za-z0-9_]", "", condition_tags$File_name)
condition_tags$Clean_Tags <- iconv(condition_tags$Clean_Tags, from = "latin1", to = "UTF-8")

significant_filtered$Clean_Tags <- condition_tags$Clean_Tags[match(significant_filtered$Condition, condition_tags$File_name)]
significant_filtered$Category <- condition_tags$Category[match(significant_filtered$Condition, condition_tags$File_name)]

significant_filtered <- significant_filtered[!is.na(significant_filtered$Category),]
significant_filtered$Count <- as.numeric(significant_filtered$Count)

significant_filtered <- significant_filtered %>%
  group_by(Clean_Tags,Condition,Category) %>%
  summarise(
    Count = sum(Count, na.rm = TRUE),
    Type = first(Type),
    Condition = first(Condition),
    Category = first(Category),
    .groups = "drop"
  )

category_colors <- c(
  "Chemical stress" = "#377EB8",
  "Carbon source utilisation" = "#4DAF4A",
  "Antibiotics" = "#E41A1C",
  "Environmental stress" = "#984EA3",
  "Antiseptic" = "#FF7F00",
  "Stress combination" = "#A65628",
  "Metal stress" = "#B8860B",
  "Base" = "#999999",
  "Envelope stress" = "#17BECF"
)


df_plot <- significant_filtered %>%
  arrange(Category) %>%
  mutate(Clean_Tags = factor(Clean_Tags, levels = unique(Clean_Tags)))

df_plot$Count <- as.numeric(df_plot$Count)

ggplot(df_plot, aes(x = Clean_Tags, y = Count, fill = Category)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = category_colors) +
  coord_flip() +
  
  labs(x = "Condition", y = "Count", fill = "Category") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 4),   # labels
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold")
  )

df_plot$indicator <- ifelse(df_plot$Count==0,"No Detected Hit","Detected Hit")
table(df_plot$indicator,df_plot$Category)

df_prop <- df_plot %>%
  count(Category, indicator) %>%
  group_by(Category) %>%
  mutate(prop = n / sum(n))

ggplot(df_prop, aes(x = Category, y = prop, fill = factor(indicator))) +
  geom_col(position = "stack", width = 0.7) +
  
  scale_fill_manual(
    values = c("No Detected Hit" = "blue", "Detected Hit" = "#E41A1C"),
    name = "GWAS Indicator"
  ) +
  
  scale_y_continuous(labels = scales::percent_format()) +
  
  theme_bw() +
  labs(x = "Category", y = "Relative frequency") +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

#Pangenome GWAS

#pangenome
dir_path <- "/Data/GWAS/panGenome"

files <- list.files(
  path = dir_path,
  pattern = ".results.csv" , 
  recursive = TRUE,
  full.names = TRUE
)

files_short <- list.files(
  path = dir_path,
  pattern = ".results.csv" , 
  recursive = TRUE
)

files_short <- gsub(".results.csv","",files_short)

significant<-c()
for(i in 1:length(files_short)){
  print(i/length(files_short))
  tmp<- read_csv(files[i],show_col_types = FALSE)
  tmp$condition<-files_short[i]
  #tmp_1<-tmp[tmp$Benjamini_H_p < 0.05 & tmp$Best_pairwise_comp_p >= 0.05 & tmp$Worst_pairwise_comp_p >= 0.05,]
  tmp_2<-tmp[tmp$Benjamini_H_p < 0.05 & tmp$Best_pairwise_comp_p < 0.05 & tmp$Worst_pairwise_comp_p >= 0.05,]
  tmp_2_unique <- unique(paste(tmp_2$Number_pos_present_in,tmp_2$Number_neg_present_in) )
  
  #tmp_3<-tmp[tmp$Benjamini_H_p < 0.05 & tmp$Best_pairwise_comp_p >= 0.05 & tmp$Worst_pairwise_comp_p < 0.05,]
  tmp_4<-tmp[tmp$Benjamini_H_p < 0.05 & tmp$Best_pairwise_comp_p < 0.05 & tmp$Worst_pairwise_comp_p < 0.05,]
  tmp_4_unique <- unique(paste(tmp_4$Number_pos_present_in,tmp_4$Number_neg_present_in) )
  
  #significant<-rbind(significant,c( "0", dim(tmp)[1],files_short[i] ))
  #significant<-rbind(significant,c( "1", dim(tmp_1)[1],files_short[i] ))
  significant<-rbind(significant,c( "Broad", length(tmp_2_unique),files_short[i] ))
  #significant<-rbind(significant,c( "3", dim(tmp_3)[1],files_short[i] ))
  significant<-rbind(significant,c( "Narrow",length(tmp_4_unique),files_short[i] ))
}

significant<- data.frame(significant)
significant$X1<- as.factor(as.character(significant$X1))
significant$X2<- as.numeric( as.character(significant$X2) )

colnames(significant) <- c("Type","Count","Condition")

significant_filtered<- significant[significant$Type == "Broad",]

significant_filtered$Condition <- gsub("[^A-Za-z0-9_]", "", significant_filtered$Condition )
significant_filtered$Condition <- gsub("_", "", significant_filtered$Condition )

condition_tags<- read_csv("/Data/condition_group.csv", locale = locale(encoding = "UTF-8"))
condition_tags$File_name <- gsub("[^A-Za-z0-9_]", "", condition_tags$File_name)
condition_tags$Clean_Tags <- iconv(condition_tags$Clean_Tags, from = "latin1", to = "UTF-8")

significant_filtered$Clean_Tags <- condition_tags$Clean_Tags[match(significant_filtered$Condition, condition_tags$File_name)]
significant_filtered$Category <- condition_tags$Category[match(significant_filtered$Condition, condition_tags$File_name)]

significant_filtered <- significant_filtered[!is.na(significant_filtered$Category),]
significant_filtered$Count <- as.numeric(significant_filtered$Count)

significant_filtered <- significant_filtered %>%
  group_by(Clean_Tags,Condition,Category) %>%
  summarise(
    Count = sum(Count, na.rm = TRUE),
    Type = first(Type),
    Condition = first(Condition),
    Category = first(Category),
    .groups = "drop"
  )

category_colors <- c(
  "Chemical stress" = "#377EB8",
  "Carbon source utilisation" = "#4DAF4A",
  "Antibiotics" = "#E41A1C",
  "Environmental stress" = "#984EA3",
  "Antiseptic" = "#FF7F00",
  "Stress combination" = "#A65628",
  "Metal stress" = "#B8860B",
  "Base" = "#999999",
  "Envelope stress" = "#17BECF"
)


df_plot <- significant_filtered %>%
  arrange(Category) %>%
  mutate(Clean_Tags = factor(Clean_Tags, levels = unique(Clean_Tags)))

df_plot$Count <- as.numeric(df_plot$Count)

ggplot(df_plot, aes(x = Clean_Tags, y = Count, fill = Category)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = category_colors) +
  coord_flip() +
  
  labs(x = "Condition", y = "Count", fill = "Category") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 4),   # labels
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold")
  )

df_plot$indicator <- ifelse(df_plot$Count==0,"No Detected Hit","Detected Hit")
table(df_plot$indicator,df_plot$Category)

df_prop <- df_plot %>%
  count(Category, indicator) %>%
  group_by(Category) %>%
  mutate(prop = n / sum(n))

ggplot(df_prop, aes(x = Category, y = prop, fill = factor(indicator))) +
  geom_col(position = "stack", width = 0.7) +
  
  scale_fill_manual(
    values = c("No Detected Hit" = "blue", "Detected Hit" = "#E41A1C"),
    name = "GWAS Indicator"
  ) +
  
  scale_y_continuous(labels = scales::percent_format()) +
  
  theme_bw() +
  labs(x = "Category", y = "Relative frequency") +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

#change 
#liyold banks offshore 
#tax advise

#no investment account 

#banking cuurrent and saving GBPs 10K worldwide 
#saving GBPs customer (saving accounts instant no fee 50K 2% fix term deposit cash deposits 6 mo 1 2 year 50K 3.25% 3.2%)
#50K usd 100K 6 month 3.3%

#mortage prduct buy UK property 1 million 30% for rent residents saudi  

#9-10 results 
#11-12 Discussion
#13-15 polish and share

#plot tree for GWAS hits 

library(ape)
library(phangorn)

pangenome_short<-read_tsv("/Data/gene_presence_absence.Rtab",show_col_types = FALSE)

colnames(pangenome_short)[which(colnames(pangenome_short)=="HA-108")]="HA.108"
colnames(pangenome_short)[which(colnames(pangenome_short)=="HA-107")]="HA.107"
colnames(pangenome_short)[which(colnames(pangenome_short)=="HA-119")]="HA.119"
colnames(pangenome_short)[which(colnames(pangenome_short)=="GR-50")]="GR.50"
colnames(pangenome_short)[which(colnames(pangenome_short)=="HA-138")]="HA.138"

tags_short <- colnames(pangenome_short)[-1]

kmer<-read_csv("/Data/total_distance_matrix.csv")
fintess_mat <- read_csv("/Data/absolute_fitness_ML.csv")
colnames(fintess_mat)[1]<-"Name"

colnames(kmer)
tags_short
fintess_mat$Name

common_names <- Reduce(intersect, list(
  colnames(kmer),
  tags_short,
  fintess_mat$Name
))


kmer_short <- kmer[match(common_names, colnames(kmer)),match(common_names, colnames(kmer))]
fintess_mat_short <- fintess_mat[match(common_names,fintess_mat$Name),]
pangenome_short_short<- pangenome_short[,c(1, match(common_names,colnames(pangenome_short)) )]


kmer_short<-as.matrix(kmer_short)*1000
kmer_short<-round(kmer_short)
kmer_short<-dist(kmer_short)

tree_tot<-nj(kmer_short)
tree_tot<-midpoint(tree_tot)
tree_tot$tip.label <- common_names

colnames(fintess_mat_short)
fintess_mat_short_choice <- fintess_mat_short[,158] %>% pull(.)
fintess_mat_short_choice <- ifelse(fintess_mat_short_choice > median(fintess_mat_short_choice),"H","L")
names(fintess_mat_short_choice) <- fintess_mat_short$Name

tree_tot
fintess_mat_short_choice

library(tidyverse)
library(ggtree)
library(treeio)
#library(treedataverse)

p <- ggtree(tree_tot) +geom_tiplab(as_ylab=TRUE, color='black', align =TRUE, linesize=.5,linetype = "dotted", hjust = 8, size=3.3)
p

feature<-as.character(fintess_mat_short_choice )
feature_df<-data.frame(list(indicator= feature)) 
row.names(feature_df)<-fintess_mat_short$Name
distinctive_colors <- c( "grey",  "black")

gheatmap(p, feature_df,
         offset = 10,
         width = 10,
         colnames = FALSE,
         color = NA) +   # 👈 removes borders
  scale_x_ggtree() +
  scale_fill_manual(breaks = unique(fintess_mat_short_choice),
                    values = distinctive_colors,
                    name = "Source")



#phylogenetics
pangenome_short<-read_tsv("/Data/gene_presence_absence.Rtab",show_col_types = FALSE)

colnames(pangenome_short)[which(colnames(pangenome_short)=="HA-108")]="HA.108"
colnames(pangenome_short)[which(colnames(pangenome_short)=="HA-107")]="HA.107"
colnames(pangenome_short)[which(colnames(pangenome_short)=="HA-119")]="HA.119"
colnames(pangenome_short)[which(colnames(pangenome_short)=="GR-50")]="GR.50"
colnames(pangenome_short)[which(colnames(pangenome_short)=="HA-138")]="HA.138"

pangenome_short_short<- pangenome_short[,c(1, match(common_names,colnames(pangenome_short)) )]
pangenome_short_short_choice<- unlist(pangenome_short_short[pangenome_short_short$Gene=="lamB",] )  
pangenome_short_short_choice<- pangenome_short_short_choice[-1]
table(pangenome_short_short_choice)

feature<-as.character(pangenome_short_short_choice )
feature_df<-data.frame(list(indicator= feature)) 
row.names(feature_df)<-names(pangenome_short_short_choice)
distinctive_colors <- c( "red",  "blue")

gheatmap(p, feature_df,
         offset = 10,
         width = 10,
         colnames = FALSE,
         color = NA) +   # 👈 removes borders
  scale_x_ggtree() +
  scale_fill_manual(breaks = unique(pangenome_short_short_choice),
                    values = distinctive_colors,
                    name = "Source")

#cluster of genes
tmp1<- unlist(pangenome_short_short[pangenome_short_short$Gene=="group_10807",] )[-1]  
tmp2<- unlist(pangenome_short_short[pangenome_short_short$Gene=="lacI~~~purR",] )[-1] 
table(tmp1, tmp2)


#
cluster_df <- read_csv("/Data/prediction_new/ML_clusters.csv")  
colnames(cluster_df)[1]<- "Name"
cluster_df_short <- cluster_df[match(common_names,cluster_df$Name),]
cluster_df_short_choice <- cluster_df_short[,30] %>% pull(.)
unique(cluster_df_short_choice)

library(randomcoloR)
feature<-as.character(cluster_df_short_choice )
feature_df<-data.frame(list(indicator= feature)) 
row.names(feature_df)<- cluster_df_short$Name
distinctive_colors <- distinctColorPalette(length(unique(cluster_df_short_choice)))

gheatmap(p, feature_df,
         offset = 10,
         width = 10,
         colnames = FALSE,
         color = NA) +   # 👈 removes borders
  scale_x_ggtree() +
  scale_fill_manual(breaks = unique(cluster_df_short_choice),
                    values = distinctive_colors,
                    name = "Source")

#show fitness 
fintess_mat <- read_csv("/Data/absolute_fitness_ML.csv")
colnames(fintess_mat)[1]<-"Name"

condition <- "LB"
isolate_id <- "DKPU071"

values <- fintess_mat[[condition]]

iso_value <- fintess_mat %>%
  dplyr::filter(Name == isolate_id) %>%
  dplyr::pull(condition)

med_value <- median(fintess_mat[[condition]], na.rm = TRUE)
# plot
ggplot(fintess_mat, aes(x = .data[[condition]])) +
  
  # continuous distribution
  geom_density(fill = "grey80", color = "black", alpha = 0.7) +
  
  # median line
  geom_vline(xintercept = med_value,
             linetype = "dotted",
             color = "black",
             size = 1) +
  
  # isolate value
  geom_vline(xintercept = iso_value,
             color = "red",
             size = 1.2) +
  
  labs(
    title = paste("Distribution:", condition),
    subtitle = paste("Isolate:", isolate_id),
    x = "Fitness",
    y = "Density"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10, angle = 0),
    axis.text.y = element_text(size = 4),   # labels
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold")
  )


#multiple isolates 
isolate_list_df  <- read_csv("/Data/gene_1.csv")
isolate_list_df$ID <- as.character( sapply(isolate_list_df$ID, function(x) strsplit(x," ")[[1]][1]))
condition <- "LB"

# 👉 your list of isolates
isolate_list <- isolate_list_df$ID

# extract isolate values
iso_df <- fintess_mat %>%
  filter(Name %in% isolate_list) %>%
  select(Name, all_of(condition)) %>%
  rename(Fitness = all_of(condition))

# median
med_value <- median(fintess_mat[[condition]], na.rm = TRUE)

# plot
ggplot(fintess_mat, aes(x = .data[[condition]])) +
  
  # distribution
  geom_density(fill = "grey80", color = "black", alpha = 0.7) +
  
  # median
  geom_vline(xintercept = med_value,
             linetype = "dotted",
             color = "black",
             size = 1) +
  
  # isolates (multiple lines)
  geom_vline(data = iso_df,
             aes(xintercept = Fitness),
             color = "red",
             alpha = 0.01,   # 👈 opacity (0 = transparent, 1 = solid)
             size = 1.2) +
  
  labs(
    title = paste("Distribution:", condition),
    x = "Fitness",
    y = "Density",
    color = "Isolate"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10, angle = 0),
    axis.text.y = element_text(size = 4),   # labels
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold")
  )


#baplots
isolate_list_df  <- read_csv("/Data/gene_1.csv")
isolate_list_df$ID <- as.character( sapply(isolate_list_df$ID, function(x) strsplit(x," ")[[1]][1]))

fintess_mat <- read_csv("/Data/absolute_fitness_ML.csv")
colnames(fintess_mat)[1]<-"Name"

condition <- "LB"
values <- fintess_mat[[condition]]

Type <- rep("Absent", length(values))
Type[match(isolate_list_df$ID, fintess_mat$Name)] <- "Present"

tmp<-fintess_mat[match(isolate_list_df$ID, fintess_mat$Name),c(1,which(colnames(fintess_mat)==condition))] 
tmp %>%
  arrange(desc(.[[2]]))

tmp %>%
  arrange(.[[2]])


df<- data.frame(list( Type=Type , values=values  ))

ggplot(df, aes(x = Type, y = values, fill = Type)) +
  geom_boxplot() +
  labs(x = "Type", y = "Fitness")+
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 13, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 13),   # labels
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold")
  )

wilcox.test(values ~ Type, data = df)


#blast for confirmation 

#for f in *.fasta; do
#makeblastdb -in "$f" -dbtype nucl
#done

#for f in *.fasta; do
#blastn \
#-query gene_query.fa \
#-db "$f" \
#-out "${f%.fasta}.blast.tsv" \
#-outfmt "6 qseqid sseqid pident length evalue bitscore" \
#-evalue 1e-5
#done


#echo -e "Isolate\tPresence" > presence_absence.tsv

#for f in *.blast.tsv; do
#isolate=${f%.blast.tsv}

#if [ -s "$f" ]; then
#echo -e "$isolate\t1" >> presence_absence.tsv
#else
#  echo -e "$isolate\t0" >> presence_absence.tsv
#fi
#done


pangenome_short<-read_tsv("/Data/gene_presence_absence.Rtab",show_col_types = FALSE)

colnames(pangenome_short)[which(colnames(pangenome_short)=="HA-108")]="HA.108"
colnames(pangenome_short)[which(colnames(pangenome_short)=="HA-107")]="HA.107"
colnames(pangenome_short)[which(colnames(pangenome_short)=="HA-119")]="HA.119"
colnames(pangenome_short)[which(colnames(pangenome_short)=="GR-50")]="GR.50"
colnames(pangenome_short)[which(colnames(pangenome_short)=="HA-138")]="HA.138"

df<-data.frame(t(pangenome_short[ pangenome_short$Gene == "group_1569",])) 
colnames(df)[1] <- "Gene"

blast <- read_tsv("/Data/group_3530_presence_absence.tsv")
blast$pan_gene <- df$Gene[match( blast$Isolate, row.names(df))]
table(blast$pan_gene, blast$Presence)


#10-11 results 
#12-13 Discussion
#14-15 final read share



#group_10680	;lacZ	hypothetical protein;beta-galactosidase;Beta-galactosidase/beta-glucuronidase
#group_10807	;lacY	MFS transporter;Lactose permease;lactose permease
#lacI~~~purR	;lacI;purR	hypothetical protein;DNA-binding transcriptional repressor LacI;DNA-binding transcriptional regulator LacI/PurR family;Lac repressor	389	140	343	590

#ybbA	ybbA;	putative ABC transporter ATP-binding protein YbbA;ABC transporter;ABC transporter ATP-binding protein


#barplots v2
pangenome_short<-read_tsv("/Data/gene_presence_absence.Rtab",show_col_types = FALSE)

colnames(pangenome_short)[which(colnames(pangenome_short)=="HA-108")]="HA.108"
colnames(pangenome_short)[which(colnames(pangenome_short)=="HA-107")]="HA.107"
colnames(pangenome_short)[which(colnames(pangenome_short)=="HA-119")]="HA.119"
colnames(pangenome_short)[which(colnames(pangenome_short)=="GR-50")]="GR.50"
colnames(pangenome_short)[which(colnames(pangenome_short)=="HA-138")]="HA.138"

isolate_list_df<-data.frame(t(pangenome_short[ pangenome_short$Gene == "group_3530",])) 
isolate_list_df$ID <- row.names(isolate_list_df)
colnames(isolate_list_df)[1] <- "Gene"

#isolate_list_df  <- read_csv("/Data/gene_1.csv")
#isolate_list_df$ID <- as.character( sapply(isolate_list_df$ID, function(x) strsplit(x," ")[[1]][1]))

fintess_mat <- read_csv("/Data/absolute_fitness_ML.csv")
colnames(fintess_mat)[1]<-"Name"

condition <- "LB"
values <- fintess_mat[[condition]]

shared <- intersect(fintess_mat$Name,isolate_list_df$ID )
fintess_mat <- fintess_mat[match(shared, fintess_mat$Name),]
isolate_list_df <- isolate_list_df[match(shared, isolate_list_df$ID),]

Type <- ifelse(isolate_list_df$Gene =="1", "Present","Absent" )

#Type <- rep("Absent", length(values))
#Type[match(isolate_list_df$ID, fintess_mat$Name)] <- "Present"

tmp<-fintess_mat[match(isolate_list_df$ID, fintess_mat$Name),c(1,which(colnames(fintess_mat)==condition))] 
tmp %>%
  arrange(desc(.[[2]]))

tmp %>%
  arrange(.[[2]])


df<- data.frame(list( Type=Type , values=values  ))

ggplot(df, aes(x = Type, y = values, fill = Type)) +
  geom_boxplot() +
  labs(x = "Type", y = "Fitness")+
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 11, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 11),   # labels
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13)
  )

wilcox.test(values ~ Type, data = df)

#importance feature over-representation


#write_csv(tot_importance, "/Data/predictive_totall_importance_1.csv")
importance<- read_csv("/Data/predictive_totall_importance_1.csv")

condition_tags<- read_csv("/Data/condition_group.csv", locale = locale(encoding = "UTF-8"),show_col_types = FALSE)
condition_tags$File_name <- gsub("[^A-Za-z0-9_]", "", condition_tags$File_name)
condition_tags$Clean_Tags <- iconv(condition_tags$Clean_Tags, from = "latin1", to = "UTF-8")


importance$category<- condition_tags$Category[match(importance$target, condition_tags$File_name)]
importance$Clean_Tags <- condition_tags$Clean_Tags[match(importance$target, condition_tags$File_name)]
importance <- importance[!is.na(importance$Clean_Tags),]



importance<- read_csv("/Data/predictive_totall_importance_1.csv")

1207/ sum(table(importance$feature_type))
importance_short <- importance[importance$rank_final==1,]
table(importance$feature)
length(importance$feature)
length( unique(importance$feature)[!grepl("Clust",unique(importance$feature))])
table(importance$feature_type)

table(importance_short$feature_type,importance_short$category)



df <- as.data.frame(table(importance_short$feature_type,
                          importance_short$category))

colnames(df) <- c("feature_type", "category", "count")
ggplot(df, aes(x = category, y = count, fill = feature_type)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  
  scale_fill_manual(values = c(
    "PAN" = "#66bd63",
    "SNP" = "#4575b4",
    "STR" = "#fdae61"
  )) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 13),   # labels
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold")
  )+
  
  labs(
    x = "Condition category",
    y = "Count"
  )





inp<-read_csv("/Data/predictive_features.csv")
colnames(inp)[1]<-"Name"

condition_tags<- read_csv("/Data/condition_group.csv", locale = locale(encoding = "UTF-8"),show_col_types = FALSE)
condition_tags$File_name <- gsub("[^A-Za-z0-9_]", "", condition_tags$File_name)
condition_tags$Clean_Tags <- iconv(condition_tags$Clean_Tags, from = "latin1", to = "UTF-8")


inp$category<- condition_tags$Category[match(inp$Name, condition_tags$File_name)]
inp$Clean_Tags <- condition_tags$Clean_Tags[match(inp$Name, condition_tags$File_name)]
inp <- inp[!is.na(inp$Clean_Tags),]
#write_csv(tot_importance, "/Data/predictive_totall_importance_1.csv")
inp$STR <- 50

unique(importance$category)
table(importance$Clean_Tags, importance$feature_type)

wide_df <- table(importance$Clean_Tags, importance$feature_type)
wide_df <- as.data.frame.matrix(wide_df)
wide_df$EXP_PAN <-  inp$Pangenome[match(row.names(wide_df), inp$Clean_Tags)]
wide_df$EXP_SNP <-  inp$SNP[match(row.names(wide_df), inp$Clean_Tags)]
wide_df$EXP_STR <-  inp$STR[match(row.names(wide_df), inp$Clean_Tags)]

wide_df$Category <-  inp$category[match(row.names(wide_df), inp$Clean_Tags)]

wide_df$enriched <- apply(wide_df, 1, function(x) {
  ratios <- c(
    PAN = x["PAN"] / x["EXP_PAN"],
    SNP = x["SNP"] / x["EXP_SNP"],
    STR = x["STR"] / x["EXP_STR"]
  )
  names(which.max(ratios))
})


wide_df <- wide_df %>%
  rowwise() %>%
  mutate(
    enriched = c("PAN", "SNP", "STR")[
      which.max(c(PAN/EXP_PAN, SNP/EXP_SNP, STR/EXP_STR))
    ]
  ) %>%
  ungroup()

table(wide_df$enriched, wide_df$Category)
