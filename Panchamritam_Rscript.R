# Packages
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(purrr)
library(ggplot2)
library(vegan)
library(ape)
library(ggrepel)
library(patchwork)
library(RColorBrewer)
library(ggpubr)
library(reshape2)
library(scales)
library(lme4)
library(lmerTest)
library(Hmisc)
library(igraph)
library(ggraph)
library(ggtext)
library(readr)
library(grid)

# File paths
data_dir <- "PATH/TO/DATA_FOLDER"

meta_file <- file.path(data_dir, "pan_meta_data.xlsx")
otu_file <- file.path(data_dir, "pan_otu_table.xlsx")
network_file <- file.path(data_dir, "network_detail_pan.xlsx")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Alpha diversity

meta_data <- read_excel(meta_file)
otu_table_raw <- read_excel(otu_file, sheet = "Sheet1")

meta_data <- meta_data %>%
  select(
    SampleID,
    Subject,
    Timepoint
  )

meta_data$Subject <- factor(
  meta_data$Subject,
  levels = paste0(
    "P",
    sort(
      as.numeric(
        gsub(
          "P",
          "",
          unique(meta_data$Subject)
        )
      )
    )
  )
)

meta_data$Timepoint <- factor(
  meta_data$Timepoint,
  levels = c(
    "0 hours",
    "6 hours",
    "12 hours"
  )
)

taxonomy_table <- otu_table_raw[,1:8]
otu_counts <- otu_table_raw[,9:38]
rownames(otu_counts) <- otu_table_raw$OTU_ID

otu_counts[] <- lapply(
  otu_counts,
  function(x)
    as.numeric(
      as.character(x)
    )
)

otu_counts[is.na(otu_counts)] <- 0

otu_transposed <- t(otu_counts)

otu_transposed <- as.data.frame(
  otu_transposed
)

shannon <- diversity(
  otu_transposed,
  index = "shannon"
)

simpson <- diversity(
  otu_transposed,
  index = "simpson"
)

richness <- specnumber(
  otu_transposed
)

evenness <- shannon /
  log(richness)

alpha_div <- data.frame(
  SampleID = rownames(
    otu_transposed
  ),
  
  Shannon = shannon,
  
  Simpson = simpson,
  
  Richness = richness,
  
  Evenness = evenness
)

merged_data <- merge(
  meta_data,
  alpha_div,
  by = "SampleID"
)

merged_data$Subject <- factor(
  merged_data$Subject,
  
  levels = paste0(
    "P",
    
    sort(
      as.numeric(
        gsub(
          "P",
          "",
          unique(
            merged_data$Subject
          )
        )
      )
    )
  )
)

merged_data$Subject_Timepoint <- paste(
  merged_data$Subject,
  merged_data$Timepoint,
  sep = "_"
)

subject_colors <- c(
  
  "#E41A1C",
  "#377EB8",
  "#4DAF4A",
  "#984EA3",
  "#FF7F00",
  "#A65628",
  "#F781BF",
  "#1B9E77",
  "#D95F02",
  "#7570B3"
)

shannon_p <- friedman.test(
  Shannon ~ Timepoint | Subject,
  data = merged_data
)$p.value

simpson_p <- friedman.test(
  Simpson ~ Timepoint | Subject,
  data = merged_data
)$p.value

richness_p <- friedman.test(
  Richness ~ Timepoint | Subject,
  data = merged_data
)$p.value

evenness_p <- friedman.test(
  Evenness ~ Timepoint | Subject,
  data = merged_data
)$p.value

shannon_label <- paste0("p = ",
  signif(shannon_p,3)
)

simpson_label <- paste0("p = ",
  signif(simpson_p,3)
)

richness_label <- paste0("p = ",
  signif(richness_p,3)
)

evenness_label <- paste0("p = ",
  signif(evenness_p,3)
)

common_theme <- theme_classic(
  base_size = 14
) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16),
    
    axis.title = element_text(
      face = "bold",
      size = 14),
    
    axis.text = element_text(
      face = "bold",
      color = "black"),
    
    legend.title = element_text(
      face = "bold"),
    
    legend.text = element_text(
      size = 11)
  )

create_alpha_plot <- function(
    y_variable,
    title_name,
    y_label,
    p_label
){
  
  ggplot(
    merged_data,
    
    aes(
      x = Timepoint,
      y = .data[[y_variable]],
      group = Subject,
      color = Subject
    )
    
  )+
    
    geom_line(
      linewidth = 1
    )+
    
    geom_point(
      size = 3
    )+
    
    geom_text(
      
      aes(label = Subject_Timepoint),
      size = 3,
      vjust = -0.8,
      show.legend = FALSE)+
    
    stat_summary(aes(group=1),
      fun=mean,
      geom="line",
      linewidth=2,
      color="black")+
    
    stat_summary(aes(group=1),
      fun=mean,
      geom="point",
      size=4,
      color="black")+
    
    scale_color_manual(
      values=subject_colors)+
    
    annotate(
      
      "text",
      x=2,
      y=max(
        merged_data[[y_variable]]
      )*1.05,
      label=p_label,
      size=5,
      fontface="bold")+
    
    common_theme+
    
    labs(title=title_name,
      x="Fermentation Time",
      y=y_label)
}

p1 <- create_alpha_plot(
  y_variable="Shannon",
  title_name="Shannon Diversity",
  y_label="Shannon Index",
  p_label=shannon_label)

p2 <- create_alpha_plot(
  y_variable="Simpson",
  title_name="Simpson Diversity",
  y_label="Simpson Index",
  p_label=simpson_label)

p3 <- create_alpha_plot(
  y_variable="Richness",
  title_name="Observed Species",
  y_label="Observed Species",
  p_label=richness_label)

p4 <- create_alpha_plot(
  y_variable="Evenness",
  title_name="Pielou Evenness",
  y_label="Pielou Evenness",
  p_label=evenness_label)

print(p1)
print(p2)
print(p3)
print(p4)

ggsave(
  file.path(output_dir, "Shannon_Diversity.png"),
  p1,
  width=8,
  height=6,
  dpi=600,
  bg="white"
)

ggsave(
  file.path(output_dir, "Simpson_Diversity.png"),
  p2,
  width=8,
  height=6,
  dpi=600,
  bg="white"
)

ggsave(
  file.path(output_dir, "Observed_Species.png"),
  p3,
  width=8,
  height=6,
  dpi=600,
  bg="white"
)

ggsave(
  file.path(output_dir, "Pielou_Evenness.png"),
  p4,
  width=8,
  height=6,
  dpi=600,
  bg="white"
)

cat("\nShannon:",shannon_label)
cat("\nSimpson:",simpson_label)
cat("\nObserved Species:",richness_label)
cat("\nPielou Evenness:",evenness_label)

# Bray-Curtis PCoA (all timepoints)

bray_dist <- vegdist(
  otu_transposed,
  method = "bray")

pcoa_result <- pcoa(
  bray_dist
)

pcoa_df <- data.frame(
  SampleID = rownames(
    otu_transposed
  ),
  
  Axis1 = pcoa_result$vectors[,1],
  Axis2 = pcoa_result$vectors[,2]
  )

beta_data <- merge(
  pcoa_df,
  meta_data,
  by = "SampleID")

beta_data$Subject <- factor(
  beta_data$Subject,
  levels = paste0(
    "P",
    sort(
      as.numeric(
        gsub("P", "",
          unique(beta_data$Subject)
        ) )) ))

beta_data$Timepoint <- factor(
  beta_data$Timepoint,
  levels = c(
    "0 hours",
    "6 hours",
    "12 hours")
)

permanova <- adonis2(
  bray_dist ~
    Timepoint,
  data = meta_data,
  permutations = 999,
  by = NULL
)

R2_value <- round(
  permanova$R2[1],
   3)

p_value <- permanova$`Pr(>F)`[1]
p_label <- paste0(
  "PERMANOVA\n",
  "R² = ",
  R2_value,
  "\n",
  "p = ",
  signif(
    p_value,
    3
  ))

axis1_percent <- round(
  pcoa_result$values$Relative_eig[1] * 100,
  1)

axis2_percent <- round(
  pcoa_result$values$Relative_eig[2] * 100,
  1)

subject_colors <- c(
  "#E41A1C",
  "#377EB8",
  "#4DAF4A",
  "#984EA3",
  "#FF7F00",
  "#A65628",
  "#F781BF",
  "#1B9E77",
  "#D95F02",
  "#7570B3"
)

beta_plot <- ggplot(
  beta_data,
  aes(
    x = Axis1,
    y = Axis2,
    color = Subject,
    shape = Timepoint
  ))+
  
  geom_point(
    size = 5
  )+
  
  geom_path(
    aes(group = Subject),
    linewidth = 1
  )+
  
  geom_text_repel(
    aes(label = Subject),
    size = 4,
    show.legend = FALSE
  )+
  
  scale_color_manual(
    values = subject_colors
  )+
  
  annotate(
    "text",
    x = max(beta_data$Axis1),
    y = max(beta_data$Axis2),
    label = p_label,
    hjust = 1,
    vjust = 1,
    fontface = "bold",
    size = 5
  )+
  
  theme_classic(base_size = 14
  )+
  
  theme(plot.title = element_text(
      
      face = "bold",
      hjust = 0.5
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    axis.text = element_text(
      face = "bold",
      color = "black"
    ),
    
    legend.title = element_text(
      face = "bold"
    ),
    
    legend.position = "right"
  )+
  
  labs(
    title = "Beta Diversity",
    x = paste0(
      "PCoA1 (", axis1_percent, "%)" ),
    
    y = paste0(
      "PCoA2 (", axis2_percent, "%)" ),
    
    color = "Subject",
    shape = "Timepoint"
  )

print(beta_plot)

ggsave(
  file.path(output_dir, "Beta_Diversity_BrayCurtis.png"),
  beta_plot,
  width = 9,
  height = 7,
  dpi = 600,
  bg = "white"
)

cat("\nPERMANOVA RESULTS\n")
cat("\nR2 =", R2_value)
cat("\np value =", signif(p_value, 3))

# Bray-Curtis PCoA (0 h vs 12 h)

otu_raw_panel <- readxl::read_excel(otu_file, sheet = "filtered_otu") |>
  as.data.frame(check.names = FALSE)
metadata_panel <- meta_data |>
  as.data.frame(check.names = FALSE)

otu_raw <- otu_raw_panel
metadata <- metadata_panel

required_metadata_columns <- c("SampleID", "Subject", "Timepoint", "Iron")
missing_metadata_columns <- setdiff(required_metadata_columns, names(metadata))
if (length(missing_metadata_columns) > 0) {
  stop("Metadata is missing required columns: ",
       paste(missing_metadata_columns, collapse = ", "))
}

required_taxonomy_columns <- c("OTU_ID", "Genus", "Species")
missing_taxonomy_columns <- setdiff(required_taxonomy_columns, names(otu_raw))
if (length(missing_taxonomy_columns) > 0) {
  stop("Filtered OTU table is missing columns: ",
       paste(missing_taxonomy_columns, collapse = ", "))
}

metadata <- metadata |>
  mutate(
    SampleID = as.character(SampleID),
    Subject = as.character(Subject),
    Timepoint = stringr::str_trim(as.character(Timepoint)),
    Iron = readr::parse_number(as.character(Iron))
  ) |>
  filter(!is.na(SampleID), SampleID != "") |>
  distinct(SampleID, .keep_all = TRUE)

sample_columns <- intersect(metadata$SampleID, names(otu_raw))
if (length(sample_columns) == 0) {
  stop("No metadata SampleID values match sample columns in filtered_otu.")
}

otu_raw[sample_columns] <- lapply(otu_raw[sample_columns], function(x) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  x[is.na(x)] <- 0
  x
})

otu_tax <- otu_raw |>
  mutate(
    Genus_clean = na_if(stringr::str_trim(as.character(Genus)), ""),
    Species_clean = na_if(stringr::str_trim(as.character(Species)), ""),
    Taxon = case_when(
      !is.na(Genus_clean) & !is.na(Species_clean) ~ paste(Genus_clean, Species_clean),
      !is.na(Genus_clean) ~ Genus_clean,
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(Taxon))

# Bray-Curtis PCoA (0 h vs 12 h)

meta_0_12 <- metadata |>
  filter(Timepoint %in% c("0 hours", "12 hours"), SampleID %in% sample_columns) |>
  mutate(
    Timepoint = factor(Timepoint, levels = c("0 hours", "12 hours"))
  )

if (nrow(meta_0_12) < 4) {
  stop("Fewer than four 0-hour/12-hour samples were found.")
}

counts_0_12 <- otu_tax |>
  select(Taxon, all_of(meta_0_12$SampleID)) |>
  group_by(Taxon) |>
  summarise(
    across(all_of(meta_0_12$SampleID), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  column_to_rownames("Taxon") |>
  as.matrix()

count_matrix_0_12 <- t(counts_0_12)
count_matrix_0_12 <- count_matrix_0_12[
  meta_0_12$SampleID, , drop = FALSE
]

count_matrix_0_12 <- count_matrix_0_12[
  , colSums(count_matrix_0_12) > 0, drop = FALSE
]

bray_0_12 <- vegan::vegdist(
  count_matrix_0_12,
  method = "bray"
)

pcoa_0_12 <- cmdscale(
  bray_0_12,
  k = 2,
  eig = TRUE,
  add = TRUE
)

pcoa_df <- as.data.frame(pcoa_0_12$points) |>
  rownames_to_column("SampleID")

names(pcoa_df)[2:3] <- c("PCoA1", "PCoA2")

pcoa_df <- pcoa_df |>
  left_join(
    meta_0_12 |>
      select(SampleID, Subject, Timepoint),
    by = "SampleID"
  )

positive_eigenvalues <- pcoa_0_12$eig[pcoa_0_12$eig > 0]
axis1_pct <- 100 * positive_eigenvalues[1] / sum(positive_eigenvalues)
axis2_pct <- 100 * positive_eigenvalues[2] / sum(positive_eigenvalues)

set.seed(123)
permanova_0_12 <- vegan::adonis2(
  bray_0_12 ~ Timepoint,
  data = meta_0_12,
  permutations = 999
)

permanova_r2 <- permanova_0_12$R2[1]
permanova_p <- permanova_0_12$`Pr(>F)`[1]

permanova_label <- paste0(
  "PERMANOVA\nR² = ", sprintf("%.3f", permanova_r2),
  "\np = ", format.pval(permanova_p, digits = 2, eps = 0.001)
)

pcoa_plot_1B <- ggplot(
  pcoa_df,
  aes(x = PCoA1, y = PCoA2, color = Timepoint)
) +
  geom_point(size = 3.5, alpha = 0.9) +
  stat_ellipse(
    aes(group = Timepoint),
    type = "t",
    linewidth = 0.8,
    level = 0.95
  ) +
  ggrepel::geom_text_repel(
    aes(label = Subject),
    size = 3,
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  annotate(
    "text",
    x = min(pcoa_df$PCoA1, na.rm = TRUE),
    y = max(pcoa_df$PCoA2, na.rm = TRUE),
    label = permanova_label,
    hjust = 0,
    vjust = 1,
    size = 3.8,
    fontface = "bold"
  ) +
  scale_color_manual(
    values = c("0 hours" = "#D98BB8", "12 hours" = "#7560A8")
  ) +
  labs(
    title = "Beta Diversity: 0 Hours vs 12 Hours",
    x = paste0("PCoA1 (", round(axis1_pct, 2), "%)"),
    y = paste0("PCoA2 (", round(axis2_pct, 2), "%)"),
    color = "Timepoint"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

print(pcoa_plot_1B)

ggsave(
  file.path(output_dir, "Bray_Curtis_PCoA_0_vs_12h.tiff"),
  pcoa_plot_1B,
  width = 7.5, height = 6, dpi = 600, compression = "lzw"
)

write.csv(
  as.data.frame(permanova_0_12),
  file.path(output_dir, "Bray_Curtis_PERMANOVA_results.csv"),
  row.names = TRUE
)

# CLR-based longitudinal LME

otu_table <- read_excel(otu_file, sheet = "filtered_otu")

meta_data <- read_excel(meta_file, sheet = "Sheet1")

meta_lme <- meta_data %>%
  filter(
    Timepoint %in% c(
      "0 hours",
      "6 hours",
      "12 hours"
    )
  ) %>%
  mutate(
    Timepoint = factor(
      Timepoint,
      levels = c(
        "0 hours",
        "6 hours",
        "12 hours"
      )
    ),
    Subject = factor(Subject)
  )

otu_clean <- otu_table %>%
  mutate(
    
    Genus = ifelse(
      is.na(Genus) | trimws(Genus) == "",
      "unclassified",
      trimws(Genus)
    ),
    
    Species = ifelse(
      is.na(Species) |
        trimws(Species) == "" |
        grepl(
          "unclassified|uncultured|unknown",
          Species,
          ignore.case = TRUE
        ),
      "unclassified",
      trimws(Species)
    ),
    
    Genus_species = tolower(
      paste(
        Genus,
        Species,
        sep = "_"
      )
    )
  )

sample_ids <- meta_lme$SampleID

otu_abundance <- otu_clean %>%
  select(
    Genus_species,
    all_of(sample_ids)
  ) %>%
  group_by(Genus_species) %>%
  summarise(
    across(
      everything(),
      ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

prevalence <- otu_abundance %>%
  select(-Genus_species) %>%
  mutate(
    across(
      everything(),
      ~ .x > 0
    )
  ) %>%
  rowSums() / ncol(
    select(otu_abundance, -Genus_species)
  )

otu_abundance <- otu_abundance %>%
  mutate(
    Prevalence = prevalence
  ) %>%
  filter(
    Prevalence >= 0.30
  ) %>%
  select(
    -Prevalence
  )

cat(
  "Number of bacteria after 30% prevalence filter:",
  nrow(otu_abundance),
  "\n"
)

count_matrix <- otu_abundance %>%
  column_to_rownames("Genus_species") %>%
  as.matrix()

relative_abundance <- sweep(
  count_matrix,
  2,
  colSums(count_matrix),
  "/"
)

pseudocount <- 1e-06

relative_abundance_pc <- relative_abundance +
  pseudocount

clr_matrix <- apply(
  relative_abundance_pc,
  2,
  function(x) {
    
    log_x <- log(x)
    
    log_x - mean(log_x)
    
  }
)

rownames(clr_matrix) <- rownames(relative_abundance)

clr_long <- as.data.frame(
  clr_matrix
) %>%
  rownames_to_column(
    "Genus_species"
  ) %>%
  pivot_longer(
    cols = -Genus_species,
    names_to = "SampleID",
    values_to = "CLR_Abundance"
  ) %>%
  left_join(
    meta_lme,
    by = "SampleID"
  )

head(clr_long)

lme_results <- clr_long %>%
  group_by(Genus_species) %>%
  group_modify(~ {
    
    dat <- .x
    
    model <- tryCatch(
      
      lmer(
        CLR_Abundance ~ Timepoint + (1 | Subject),
        data = dat,
        REML = FALSE
      ),
      
      error = function(e) NULL
    )
    
    if (is.null(model)) {
      
      return(
        data.frame(
          Effect = NA,
          Estimate = NA,
          SE = NA,
          p_value = NA
        )
      )
    }
    
    coef_table <- summary(model)$coefficients
    
    result <- data.frame(
      
      Effect = rownames(coef_table),
      
      Estimate = coef_table[, "Estimate"],
      
      SE = coef_table[, "Std. Error"],
      
      p_value = coef_table[, "Pr(>|t|)"],
      
      row.names = NULL
    )
    
    result %>%
      filter(
        Effect %in% c(
          "Timepoint6 hours",
          "Timepoint12 hours"
        )
      )
    
  }) %>%
  ungroup()

lme_results <- lme_results %>%
  mutate(
    FDR = p.adjust(
      p_value,
      method = "BH"
    )
  )

significant_bacteria <- lme_results %>%
  filter(
    !is.na(p_value),
    p_value < 0.05
  ) %>%
  arrange(
    p_value
  )

significant_bacteria

plot_data <- significant_bacteria %>%
  filter(
    !is.na(Estimate),
    !is.na(p_value)
  ) %>%
  mutate(
    
    Timepoint = case_when(
      Effect == "Timepoint6 hours"  ~ "0 hours vs 6 hours",
      Effect == "Timepoint12 hours" ~ "0 hours vs 12 hours"
    ),
    
    Bacteria = gsub(
      "_",
      " ",
      Genus_species
    ),
    
    Bacteria = paste0(
      Bacteria
    ),
    
    Direction = ifelse(
      Estimate > 0,
      "Increased",
      "Decreased"
    )
  )

ggplot(
  plot_data,
  aes(
    x = Estimate,
    y = reorder(Bacteria, Estimate),
    color = Direction
  )
) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  geom_point(
    size = 4
  ) +
  
  geom_errorbarh(
    aes(
      xmin = Estimate - SE,
      xmax = Estimate + SE
    ),
    height = 0.25,
    linewidth = 0.8
  ) +
  
  facet_wrap(
    ~ Timepoint,
    ncol = 1
  ) +
  
  scale_color_manual(
    values = c(
      "Increased" = "#4575B4",
      "Decreased" = "#D73027"
    )
  ) +
  
  labs(
    title = "Bacterial Taxa Associated with Timepoint",
    subtitle = "Linear mixed-effects model | 0 hours as reference",
    x = "LME coefficient estimate ± SE (CLR abundance)",
    y = NULL,
    color = "Direction"
  ) +
  
  theme_classic() +
  
  theme(
    
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 15
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    
    axis.text.y = element_text(
      face = "italic",
      size = 10
    ),
    
    axis.title.x = element_text(
      face = "bold"
    ),
    
    legend.position = "top"
  )

# P8 mean relative-abundance difference

p8_meta <- metadata |>
  filter(
    Subject == "P8",
    Timepoint %in% c("0 hours", "6 hours"),
    SampleID %in% sample_columns
  ) |>
  mutate(
    Timepoint = factor(Timepoint, levels = c("0 hours", "6 hours"))
  )

if (nrow(p8_meta) == 0) {
  stop("No P8 samples at 0 or 6 hours were found.")
}

p8_species_for_plot <- c(
  "Blautia obeum",
  "Bifidobacterium longum",
  "Catenibacterium sp. co_0103",
  "Subdoligranulum variabile",
  "Mediterraneibacter faecis",
  "Mediterraneibacter lactaris",
  "Bifidobacterium adolescentis",
  "Ligilactobacillus ruminis",
  "Xanthomonas albilineans",
  "Streptococcus thermophilus"
)

p8_counts_long <- otu_tax |>
  select(Taxon, all_of(p8_meta$SampleID)) |>
  pivot_longer(
    cols = all_of(p8_meta$SampleID),
    names_to = "SampleID",
    values_to = "Count"
  ) |>
  group_by(SampleID) |>
  mutate(
    Total_reads = sum(Count, na.rm = TRUE),
    Relative_Abundance = ifelse(
      Total_reads > 0, Count / Total_reads, NA_real_
    )
  ) |>
  ungroup() |>
  left_join(
    p8_meta |>
      select(SampleID, Subject, Timepoint),
    by = "SampleID"
  ) |>
  group_by(Taxon, SampleID, Subject, Timepoint) |>
  summarise(
    Relative_Abundance = sum(Relative_Abundance, na.rm = TRUE),
    .groups = "drop"
  )

p8_abundance_complete <- tidyr::expand_grid(
  Taxon = p8_species_for_plot,
  SampleID = p8_meta$SampleID
) |>
  left_join(
    p8_counts_long,
    by = c("Taxon", "SampleID")
  ) |>
  left_join(
    p8_meta |>
      select(SampleID, Subject, Timepoint),
    by = "SampleID",
    suffix = c("", ".meta")
  ) |>
  mutate(
    Timepoint = coalesce(Timepoint, Timepoint.meta),
    Relative_Abundance = replace_na(Relative_Abundance, 0)
  ) |>
  select(Taxon, SampleID, Subject, Timepoint, Relative_Abundance)

p8_mean_abundance <- p8_abundance_complete |>
  group_by(Taxon, Timepoint) |>
  summarise(
    Mean_Relative_Abundance = mean(Relative_Abundance, na.rm = TRUE),
    .groups = "drop"
  ) |>
  pivot_wider(
    names_from = Timepoint,
    values_from = Mean_Relative_Abundance,
    values_fill = 0
  ) |>
  mutate(
    Difference_6h_minus_0h = `6 hours` - `0 hours`
  ) |>
  arrange(desc(Difference_6h_minus_0h))

p8_mean_abundance <- p8_mean_abundance |>
  mutate(
    Taxon = factor(Taxon, levels = Taxon)
  )

p8_plot_2B <- ggplot(
  p8_mean_abundance,
  aes(
    x = Taxon,
    y = Difference_6h_minus_0h,
    fill = Difference_6h_minus_0h > 0
  )
) +
  geom_col(
    width = 0.8,
    color = "black",
    linewidth = 0.4
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  scale_fill_manual(
    values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"),
    guide = "none"
  ) +
  scale_x_discrete(
    labels = function(x) gsub(" ", "_", x)
  ) +
  labs(
    title = "Difference in Mean Relative Abundance (6 h - 0 h)",
    x = "Species",
    y = "Difference in Relative Abundance"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1,
      face = "bold", color = "black"
    )
  )

print(p8_plot_2B)

ggsave(
  file.path(output_dir, "P8_mean_relative_abundance_difference.tiff"),
  p8_plot_2B,
  width = 10, height = 6.5, dpi = 600, compression = "lzw"
)

write.csv(
  p8_mean_abundance,
  file.path(output_dir, "P8_mean_abundance_difference.csv"),
  row.names = FALSE
)

# SCFA-associated bacteria
meta_data <- readxl::read_excel(meta_file, sheet = "Sheet1")
otu_table <- readxl::read_excel(otu_file, sheet = "filtered_otu")

meta_scfa <- meta_data %>%
  filter(
    Timepoint %in% c(
      "0 hours",
      "6 hours",
      "12 hours"
    )
  ) %>%
  mutate(
    Timepoint = factor(
      Timepoint,
      levels = c(
        "0 hours",
        "6 hours",
        "12 hours"
      )
    )
  )

sample_ids <- meta_scfa$SampleID

otu_clean <- otu_table %>%
  mutate(
    
    Genus = ifelse(
      is.na(Genus) |
        trimws(Genus) == "",
      "unclassified",
      trimws(Genus)
    ),
    
    Species = ifelse(
      is.na(Species) |
        trimws(Species) == "",
      "unclassified",
      trimws(Species)
    ),
    
    Species = ifelse(
      grepl(
        "unclassified|uncultured|unknown",
        Species,
        ignore.case = TRUE
      ),
      "unclassified",
      Species
    ),
    
    Genus_species = tolower(
      paste(
        Genus,
        Species,
        sep = "_"
      )
    ),
    
    Genus_clean = tolower(
      trimws(Genus)
    )
  )

acetate_specific <- c(
  "akkermansia_muciniphila",
  "clostridium_ramosum"
)

acetate_genera <- c(
  "bacteroides",
  "bifidobacterium",
  "clostridium",
  "desulfovibrio",
  "lactobacillus",
  "prevotella",
  "ruminococcus"
)

propionate_specific <- c(
  "alistipes_putredinis",
  "blautia_obeum",
  "clostridium_ramosum",
  "coprococcus_catus",
  "eubacterium_hallii",
  "megasphaera_elsdenii",
  "phascolarcotobacterium_succinatutens",
  "roseburia_inulivorans",
  "ruminococcus_obeum"
)

propionate_genera <- c(
  "bacteroides",
  "dialister",
  "salmonella",
  "veillonella"
)

butyrate_specific <- c(
  "anaerostipes_caccae",
  "clostridium_leptum",
  "coprococcus_catus",
  "coprococcus_comes",
  "coprococcus_eutacus",
  "eubacterium_biforme",
  "eubacterium_hallii",
  "eubacterium_rectale",
  "faecalibacterium_prausnitzii",
  "roseburia_inulivorans",
  "subdoligranulum_variabile"
)

butyrate_genera <- c(
  "anaerostipes",
  "holdemanella",
  "roseburia"
)

otu_scfa <- otu_clean %>%
  mutate(

    Acetate = ifelse(
      Genus_species %in% acetate_specific |
        Genus_clean %in% acetate_genera,
      1,
      0
    ),

    Propionate = ifelse(
      Genus_species %in% propionate_specific |
        Genus_clean %in% propionate_genera,
      1,
      0
    ),

    Butyrate = ifelse(
      Genus_species %in% butyrate_specific |
        Genus_clean %in% butyrate_genera,
      1,
      0
    )
  )

scfa_bacteria_present <- otu_scfa %>%
  filter(
    Acetate == 1 |
      Propionate == 1 |
      Butyrate == 1
  ) %>%
  select(
    Genus_species,
    Genus_clean,
    Acetate,
    Propionate,
    Butyrate
  ) %>%
  distinct()

cat(
  "\nNumber of SCFA-associated taxa found:",
  n_distinct(
    scfa_bacteria_present$Genus_species
  ),
  "\n"
)

print(scfa_bacteria_present)

scfa_abundance <- otu_scfa %>%
  filter(
    Acetate == 1 |
      Propionate == 1 |
      Butyrate == 1
  ) %>%
  select(
    Genus_species,
    all_of(sample_ids),
    Acetate,
    Propionate,
    Butyrate
  ) %>%
  group_by(
    Genus_species
  ) %>%
  summarise(
    
    across(
      all_of(sample_ids),
      ~ sum(.x, na.rm = TRUE)
    ),
    
    Acetate = max(
      Acetate,
      na.rm = TRUE
    ),
    
    Propionate = max(
      Propionate,
      na.rm = TRUE
    ),
    
    Butyrate = max(
      Butyrate,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

print(
  scfa_abundance %>%
    select(
      Genus_species,
      Acetate,
      Propionate,
      Butyrate
    )
)

abundance_matrix <- scfa_abundance %>%
  select(
    Genus_species,
    all_of(sample_ids)
  ) %>%
  column_to_rownames(
    "Genus_species"
  ) %>%
  as.matrix()

sample_totals <- otu_table %>%
  select(
    all_of(sample_ids)
  ) %>%
  summarise(
    across(
      everything(),
      ~ sum(.x, na.rm = TRUE)
    )
  )

sample_totals <- as.numeric(
  sample_totals[1, ]
)

names(sample_totals) <- sample_ids

relative_abundance <- sweep(
  abundance_matrix,
  2,
  sample_totals,
  "/"
) * 100

rel_abund_long <- as.data.frame(
  relative_abundance
) %>%
  rownames_to_column(
    "Genus_species"
  ) %>%
  pivot_longer(
    cols = -Genus_species,
    names_to = "SampleID",
    values_to = "Relative_Abundance"
  ) %>%
  left_join(
    scfa_abundance %>%
      select(
        Genus_species,
        Acetate,
        Propionate,
        Butyrate
      ),
    by = "Genus_species"
  ) %>%
  left_join(
    meta_scfa,
    by = "SampleID"
  )

scfa_long <- rel_abund_long %>%
  pivot_longer(
    cols = c(
      Acetate,
      Propionate,
      Butyrate
    ),
    names_to = "SCFA",
    values_to = "Producer"
  ) %>%
  filter(
    Producer == 1
  ) %>%
  group_by(
    SampleID,
    Subject,
    Timepoint,
    SCFA
  ) %>%
  summarise(
    Total_Relative_Abundance =
      sum(
        Relative_Abundance,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

scfa_summary <- scfa_long %>%
  group_by(
    SCFA,
    Timepoint
  ) %>%
  summarise(
    
    Mean = mean(
      Total_Relative_Abundance,
      na.rm = TRUE
    ),
    
    SD = sd(
      Total_Relative_Abundance,
      na.rm = TRUE
    ),
    
    N = n(),
    
    .groups = "drop"
  ) %>%
  mutate(
    Timepoint = factor(
      Timepoint,
      levels = c(
        "0 hours",
        "6 hours",
        "12 hours"
      )
    )
  )

wilcox_scfa <- scfa_long %>%
  group_by(SCFA) %>%
  group_modify(~ {
    
    x_0 <- .x %>%
      filter(
        Timepoint == "0 hours"
      ) %>%
      pull(
        Total_Relative_Abundance
      )
    
    x_6 <- .x %>%
      filter(
        Timepoint == "6 hours"
      ) %>%
      pull(
        Total_Relative_Abundance
      )
    
    x_12 <- .x %>%
      filter(
        Timepoint == "12 hours"
      ) %>%
      pull(
        Total_Relative_Abundance
      )
    
    p_0_6 <- wilcox.test(
      x_0,
      x_6,
      paired = FALSE,
      exact = FALSE
    )$p.value
    
    p_6_12 <- wilcox.test(
      x_6,
      x_12,
      paired = FALSE,
      exact = FALSE
    )$p.value
    
    data.frame(
      p_0_6 = p_0_6,
      p_6_12 = p_6_12
    )
    
  }) %>%
  ungroup()

final_scfa_table <- scfa_summary %>%
  mutate(
    Mean_SD = paste0(
      round(Mean, 2),
      " ± ",
      round(SD, 2)
    )
  ) %>%
  select(
    SCFA,
    Timepoint,
    Mean_SD
  ) %>%
  pivot_wider(
    names_from = Timepoint,
    values_from = Mean_SD
  ) %>%
  left_join(
    wilcox_scfa,
    by = "SCFA"
  ) %>%
  mutate(
    p_0_6 = sprintf(
      "%.3f",
      p_0_6
    ),
    p_6_12 = sprintf(
      "%.3f",
      p_6_12
    )
  )

print(final_scfa_table)

ggplot(
  scfa_summary,
  aes(
    x = Timepoint,
    y = Mean,
    group = SCFA,
    color = SCFA
  )
) +
  
  geom_line(
    linewidth = 1.5
  ) +
  
  scale_color_manual(
    values = c(
      "Acetate" = "steelblue",
      "Propionate" = "hotpink2",
      "Butyrate" = "#757"
    )
  ) +
  
  scale_x_discrete(
    limits = c(
      "0 hours",
      "6 hours",
      "12 hours"
    )
  ) +
  
  labs(
    title = "Changes in SCFA-Producing Bacteria Over Time",
    x = "Timepoint",
    y = "Mean relative abundance (%)",
    color = "SCFA"
  ) +
  
  theme_classic() +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 15
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    legend.title = element_text(
      face = "bold"
    ),
    
    legend.position = "right"
  )

scfa_composition <- rel_abund_long %>%
  
  pivot_longer(
    cols = c(
      Acetate,
      Propionate,
      Butyrate
    ),
    names_to = "SCFA",
    values_to = "Producer"
  ) %>%
  
  filter(
    Producer == 1
  ) %>%
  
  group_by(
    SCFA,
    Genus_species,
    Timepoint
  ) %>%
  summarise(
    Mean_Abundance = mean(
      Relative_Abundance,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  
  mutate(
    Timepoint = factor(
      Timepoint,
      levels = c(
        "0 hours",
        "6 hours",
        "12 hours"
      )
    )
  )

acetate_plot <- scfa_composition %>%
  filter(
    SCFA == "Acetate"
  ) %>%
  
  ggplot(
    aes(
      x = Timepoint,
      y = Mean_Abundance,
      fill = Genus_species
    )
  ) +
  
  geom_col(
    position = "stack",
    width = 0.7
  ) +
  
  labs(
    title = "Composition of Acetate-Producing Bacteria",
    x = "Timepoint",
    y = "Mean relative abundance (%)",
    fill = "Bacterium"
  ) +
  
  theme_classic() +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 15
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    axis.text = element_text(
      size = 10
    ),
    
    legend.title = element_text(
      face = "bold"
    ),
    
    legend.text = element_text(
      face = "italic"
    )
  )

acetate_plot

butyrate_plot <- scfa_composition %>%
  filter(
    SCFA == "Butyrate"
  ) %>%
  
  ggplot(
    aes(
      x = Timepoint,
      y = Mean_Abundance,
      fill = Genus_species
    )
  ) +
  
  geom_col(
    position = "stack",
    width = 0.7
  ) +
  
  labs(
    title = "Composition of Butyrate-Producing Bacteria",
    x = "Timepoint",
    y = "Mean relative abundance (%)",
    fill = "Bacterium"
  ) +
  
  theme_classic() +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 15
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    axis.text = element_text(
      size = 10
    ),
    
    legend.title = element_text(
      face = "bold"
    ),
    
    legend.text = element_text(
      face = "italic"
    )
  )

butyrate_plot

propionate_plot <- scfa_composition %>%
  filter(
    SCFA == "Propionate"
  ) %>%
  
  ggplot(
    aes(
      x = Timepoint,
      y = Mean_Abundance,
      fill = Genus_species
    )
  ) +
  
  geom_col(
    position = "stack",
    width = 0.7
  ) +
  
  labs(
    title = "Composition of Propionate-Producing Bacteria",
    x = "Timepoint",
    y = "Mean relative abundance (%)",
    fill = "Bacterium"
  ) +
  
  theme_classic() +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 15
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    axis.text = element_text(
      size = 10
    ),
    
    legend.title = element_text(
      face = "bold"
    ),
    
    legend.text = element_text(
      face = "italic"
    )
  )

propionate_plot

# Vitamin B-associated bacteria

vitamin_bacteria <- tribble(
  
  ~Genus_species,                         ~Vitamin_B1, ~Vitamin_B2, ~Vitamin_B3, ~Vitamin_B5, ~Vitamin_B6, ~Vitamin_B7, ~Vitamin_B9, ~Vitamin_B12,
  
  "prevotella_unclassified",                 1, 0, 0, 0, 0, 0, 0, 0,
  "desulfovibrio_unclassified",              1, 0, 0, 0, 0, 0, 0, 0,
  "bacteroides_fragilis",                    1, 1, 1, 1, 1, 1, 1, 1,
  "prevotella_copri",                        1, 1, 1, 1, 1, 1, 1, 1,
  "clostridium_difficile",                   1, 1, 1, 0, 0, 0, 1, 1,
  "ruminococcus_lactaris",                   1, 0, 1, 1, 0, 0, 0, 1,
  "lactobacillus_casei",                     1, 0, 0, 0, 0, 0, 0, 0,
  "bifidobacterium_bifidum",                 1, 0, 0, 0, 0, 0, 1, 0,
  "bifidobacterium_infantis",                1, 0, 1, 0, 0, 0, 0, 1,
  "fusobacterium_varium",                    1, 0, 1, 0, 0, 1, 1, 1,
  
  "lactobacillus_plantarum",                 0, 1, 0, 0, 0, 0, 1, 1,
  "lactobacillus_fermentum",                 0, 1, 0, 0, 0, 0, 0, 0,
  "lactobacillus_acidophilus",               0, 1, 0, 0, 0, 0, 0, 0,
  "mycobacterium_tuberculosis",              0, 1, 0, 0, 0, 0, 0, 0,
  "salmonella_typhimurium",                  0, 1, 0, 1, 0, 0, 0, 0,
  
  "helicobacter_pylori",                     0, 0, 1, 1, 1, 0, 0, 0,
  
  "escherichia_coli",                        0, 0, 0, 1, 0, 0, 0, 0,
  "corynebacterium_glutamicum",              0, 0, 0, 1, 0, 0, 0, 0,
  "ruminococcus_torques",                    0, 0, 0, 1, 0, 0, 0, 0,
  "salmonella_enterica",                     0, 0, 0, 1, 0, 0, 1, 0,
  
  "bifidobacterium_longum",                  0, 0, 0, 0, 1, 0, 1, 1,
  "collinsella_aerofaciens",                 0, 0, 0, 0, 1, 0, 0, 0,
  
  "campylobacter_coli",                      0, 0, 0, 0, 0, 1, 0, 0,
  "lactobacillus_helveticus",                0, 0, 0, 0, 0, 1, 1, 0,
  
  "lactobacillus_reuteri",                   0, 0, 0, 0, 0, 0, 1, 1,
  "lactobacillus_delbrueckii",               0, 0, 0, 0, 0, 0, 1, 0,
  "streptococcus_thermophilus",              0, 0, 0, 0, 0, 0, 1, 0,
  "bifidobacterium_breve",                   0, 0, 0, 0, 0, 0, 1, 0,
  "bifidobacterium_adolescentis",            0, 0, 0, 0, 0, 0, 1, 0,
  "lactobacillus_sakei",                     0, 0, 0, 0, 0, 0, 1, 0,
  "bifidobacterium_pseudocatenulatum",       0, 0, 0, 0, 0, 0, 1, 0,
  
  "pseudomonas_denitrificans",               0, 0, 0, 0, 0, 0, 0, 1,
  "bacillus_megaterium",                     0, 0, 0, 0, 0, 0, 0, 1,
  "propionibacterium_freudenreichii",        0, 0, 0, 0, 0, 0, 0, 1,
  "faecalibacterium_prausnitzii",            0, 0, 0, 0, 0, 0, 0, 1,
  "bifidobacterium_animalis",                0, 0, 0, 0, 0, 0, 0, 1,
  "lactobacillus_coryniformis",              0, 0, 0, 0, 0, 0, 0, 1
)

vitamin_bacteria

nrow(vitamin_bacteria)

vitamin_bacteria %>%
  count(Genus_species) %>%
  filter(n > 1)

otu_clean <- otu_table %>%
  mutate(
    
    Genus = ifelse(
      is.na(Genus) | trimws(Genus) == "",
      "unclassified",
      trimws(Genus)
    ),
    
    Species = ifelse(
      is.na(Species) | trimws(Species) == "",
      "unclassified",
      trimws(Species)
    ),
    
    Species = ifelse(
      grepl(
        "unclassified|uncultured|unknown",
        Species,
        ignore.case = TRUE
      ),
      "unclassified",
      Species
    ),
    
    Genus_species = tolower(
      paste(
        Genus,
        Species,
        sep = "_"
      )
    )
  )

meta_vitamin <- meta_data %>%
  filter(
    Timepoint %in% c(
      "0 hours",
      "6 hours",
      "12 hours"
    )
  ) %>%
  mutate(
    Timepoint = factor(
      Timepoint,
      levels = c(
        "0 hours",
        "6 hours",
        "12 hours"
      )
    )
  )

sample_ids <- meta_vitamin$SampleID

vitamin_bacteria_present <- vitamin_bacteria %>%
  filter(
    Genus_species %in% otu_clean$Genus_species
  )

print(vitamin_bacteria_present)

abundance_table <- otu_clean %>%
  filter(
    Genus_species %in%
      vitamin_bacteria_present$Genus_species
  ) %>%
  select(
    Genus_species,
    all_of(sample_ids)
  ) %>%
  group_by(
    Genus_species
  ) %>%
  summarise(
    across(
      all_of(sample_ids),
      ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

abundance_matrix <- abundance_table %>%
  column_to_rownames(
    "Genus_species"
  ) %>%
  as.matrix()

sample_totals <- otu_table %>%
  select(
    all_of(sample_ids)
  ) %>%
  summarise(
    across(
      everything(),
      ~ sum(.x, na.rm = TRUE)
    )
  )

sample_totals <- as.numeric(
  sample_totals[1, ]
)

names(sample_totals) <- sample_ids

relative_abundance <- sweep(
  abundance_matrix,
  2,
  sample_totals,
  "/"
) * 100

rel_abund_long <- as.data.frame(
  relative_abundance
) %>%
  rownames_to_column(
    "Genus_species"
  ) %>%
  pivot_longer(
    cols = -Genus_species,
    names_to = "SampleID",
    values_to = "Relative_Abundance"
  ) %>%
  left_join(
    vitamin_bacteria_present,
    by = "Genus_species"
  ) %>%
  left_join(
    meta_vitamin,
    by = "SampleID"
  )

vitamin_long <- rel_abund_long %>%
  pivot_longer(
    cols = starts_with("Vitamin_B"),
    names_to = "Vitamin",
    values_to = "Producer"
  ) %>%
  filter(
    Producer == 1
  ) %>%
  group_by(
    SampleID,
    Subject,
    Timepoint,
    Vitamin
  ) %>%
  summarise(
    Total_Relative_Abundance =
      sum(
        Relative_Abundance,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

vitamin_summary <- vitamin_long %>%
  group_by(
    Vitamin,
    Timepoint
  ) %>%
  summarise(
    Mean = mean(
      Total_Relative_Abundance,
      na.rm = TRUE
    ),
    SD = sd(
      Total_Relative_Abundance,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

vitamin_summary

ggplot(
  vitamin_summary,
  aes(
    x = Timepoint,
    y = Mean,
    group = Vitamin,
    color = Vitamin
  )
) +
  
  geom_line(
    linewidth = 1.2
  ) +
  
  geom_point(
    size = 3
  ) +
  
  labs(
    title = "Changes in Vitamin B-Producing Bacteria Over Time",
    subtitle = "Total relative abundance of identified vitamin-associated bacteria",
    x = "Timepoint",
    y = "Mean relative abundance (%)",
    color = "Vitamin"
  ) +
  
  theme_classic() +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 15
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    axis.title = element_text(
      face = "bold"
    ),
    axis.text = element_text(
      size = 11
    ),
    legend.title = element_text(
      face = "bold"
    )
  )

# Bray-Curtis PCoA (P8 removed)

head(meta_data)
meta_data <- read_excel(
  meta_file)
otu_table_raw <- read_excel(
  otu_file,
  sheet = "Sheet1")
meta_data <- meta_data[, c(
  "SampleID",
  "Subject",
  "Timepoint")]
meta_data <- meta_data[meta_data$Subject != "P8",]
meta_data$Subject <- factor(
  meta_data$Subject,
  levels = paste0(
    "P",
    sort(
      as.numeric(
        gsub(
          "P",
          "",
          unique(meta_data$Subject))))))

meta_data$Timepoint <- factor(
  meta_data$Timepoint,
  levels = c(
    "0 hours",
    "6 hours",
    "12 hours"))

otu_counts <- otu_table_raw[, 9:38]
rownames(otu_counts) <- otu_table_raw$OTU_ID
otu_counts[] <- lapply(
  otu_counts,
  function(x) as.numeric(as.character(x)))
otu_counts[is.na(otu_counts)] <- 0
otu_transposed <- t(otu_counts)
otu_transposed <- as.data.frame(otu_transposed)
otu_transposed <- otu_transposed[
  rownames(otu_transposed) %in% meta_data$SampleID,
  ,
  drop = FALSE]
otu_transposed <- otu_transposed[
  match(
    meta_data$SampleID,
    rownames(otu_transposed)),
  ,
  drop = FALSE]
bray_dist <- vegdist(
  otu_transposed,
  method = "bray")

pcoa_result <- cmdscale(
  bray_dist,
  k = 2,
  eig = TRUE,
  add = TRUE
)

pcoa_df <- data.frame(
  SampleID = rownames(otu_transposed),
  Axis1 = pcoa_result$points[, 1],
  Axis2 = pcoa_result$points[, 2])

beta_data <- merge(
  pcoa_df,
  meta_data,
  by = "SampleID")
beta_data$Subject <- factor(
  beta_data$Subject,
  levels = paste0(
    "P",
    sort(
      as.numeric(
        gsub(
          "P",
          "",
          unique(beta_data$Subject))))))
beta_data$Timepoint <- factor(
  beta_data$Timepoint,
  levels = c(
    "0 hours",
    "6 hours",
    "12 hours"))
permanova <- adonis2(
  bray_dist ~ Timepoint,
  data = meta_data,
  permutations = 999)
R2_value <- round(
  permanova$R2[1], 3)
p_value <- permanova$`Pr(>F)`[1]
p_label <- paste0(
  "PERMANOVA\n",
  "R² = ",
  R2_value,
  "\n",
  "p = ",
  signif(p_value, 3))

eig <- pcoa_result$eig
positive_eig <- eig[eig > 0]
axis1_percent <- round(
  positive_eig[1] /
    sum(positive_eig) * 100,  1)

axis2_percent <- round(
  positive_eig[2] /
    sum(positive_eig) * 100, 1)

subject_colors <- c(
  "#E41A1C",
  "#377EB8",
  "#4DAF4A",
  "#984EA3",
  "#FF7F00",
  "#A65628",
  "#F781BF",
  "#1B9E77",
  "#D95F02",
  "#7570B3")

beta_plot <- ggplot(
  beta_data,
  aes(
    x = Axis1,
    y = Axis2,
    color = Subject,
    shape = Timepoint)) +
  
  geom_point(size = 5) +
  
  geom_path(
    aes(group = Subject),
    linewidth = 1) +
  
  geom_text_repel(
    aes(label = Subject),
    size = 4,
    show.legend = FALSE) +
  
  scale_color_manual(values = subject_colors) +
  
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)) +
  
  annotate(
    "text",
    x = mean(range(beta_data$Axis1)),
    y = max(beta_data$Axis2) - 0.01,
    label = p_label,
    hjust = 0.5,
    vjust = 1,
    fontface = "bold",
    size = 5) +
  
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(
      face = "bold",
      color = "black"),
    
    legend.title = element_text(face = "bold"),
    legend.position = "right") +
  
  labs(title = "Bray-Curtis PCoA",
    
    x = paste0(
      "PCoA1 (",
      axis1_percent,
      "%)"),
    
    y = paste0(
      "PCoA2 (",
      axis2_percent,
      "%)"),
    shape = "Timepoint",
    color = "Subject")

print(beta_plot)
ggsave(
  file.path(output_dir, "Beta_remove_P8.png"),
  beta_plot,
  width = 9,
  height = 7,
  dpi = 600,
  bg = "white")
cat("\n====================================")
cat("\nBRAY-CURTIS PCoA — P8 REMOVED")
cat("\n====================================")
cat("\n\nSamples included:", nrow(beta_data))
cat("\nSubjects included:", length(unique(beta_data$Subject)))
cat("\n\nPCoA1:", axis1_percent, "%")
cat("\nPCoA2:", axis2_percent, "%")
cat("\n\nPERMANOVA RESULTS")
cat("\n------------------------------")
print(permanova)
cat("\nR² =", R2_value)
cat("\np value =", signif(p_value, 3))
cat("\n\nSamples per subject/timepoint:")
print(
  table(beta_data$Subject, beta_data$Timepoint))

# Top 20 species relative abundance
meta_data <- readxl::read_excel(meta_file, sheet = "Sheet1")
otu_table_raw <- readxl::read_excel(otu_file, sheet = "Sheet1")
taxonomy_table <- otu_table_raw[, 1:8]
otu_counts <- otu_table_raw[, 9:38]
rownames(otu_counts) <- otu_table_raw$OTU_ID
otu_counts[] <- lapply(otu_counts, function(x) as.numeric(as.character(x)))
otu_counts[is.na(otu_counts)] <- 0

relative_abundance <- sweep(
  otu_counts, 2,
  colSums(otu_counts),
  FUN = "/"
)

taxonomy_table$Genus_species <- paste(
  taxonomy_table$Genus,
  taxonomy_table$Species,
  sep = "_"
)

relative_abundance$Genus_species <- taxonomy_table$Genus_species

species_mean <- relative_abundance %>%
  group_by(Genus_species) %>%
  summarise(
    across(
      where(is.numeric),
      mean,
      na.rm = TRUE
    )
  )

species_mean$total_abundance <- rowSums(
  species_mean[,2:ncol(species_mean)],
  na.rm = TRUE
)

top20_species <- species_mean %>%
  arrange(
    desc(total_abundance)
  ) %>%
  slice(1:20)

top20_long <- melt(
  top20_species,
  id.vars = "Genus_species",
  variable.name = "SampleID",
  value.name = "RelativeAbundance"
)

top20_long <- top20_long %>%
  
  filter(SampleID != "total_abundance")

top20_long <- merge(
  top20_long,
  meta_data,
  by = "SampleID"
)

top20_long$Timepoint <- factor(
  top20_long$Timepoint,
  levels = c(
    "0 hours",
    "6 hours",
    "12 hours"
  )
)

top20_long$Subject <- factor(
  top20_long$Subject,
  levels = paste0(
    "P",
    sort(
      as.numeric(
        gsub(
          "P",
          "",
          unique(top20_long$Subject)
        )
      )
    )
  )
)

distinct_colors <- c(
    "#800000",
    "#FC1490",
    "#4363D8",
    "#FF7043",
    "#3CB44B",
    "#FABEBE",
    "#911EB4",
    "#46F0F0",
    "#FF2000",
    "#E6BEFF",
    "#004D40",
    "#377EB8",
    "#B8860B",
    "#1B9E77",
    "#FFD8B1",
    "#FFF700",
    "#FF00FF",
    "#7FFF00",
    "#E6194B",                             
    "#000075"
)

species_plot <- ggplot(
  top20_long,
  aes(
    x = Subject,
    y = RelativeAbundance,
    fill = Genus_species
  )
) +
  
  geom_bar(
    stat = "identity",
    position = "fill",
    width = 0.9
  ) +
  
  facet_grid(
    ~ Timepoint,
    scales = "free_x",
    space = "free_x"
  ) +
  
  scale_fill_manual(
    values = distinct_colors
  ) +
  
  scale_y_continuous(
    labels = percent_format()
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 11,
      face = "bold",
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 11,
      face = "bold",
      color = "black"
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 14
    ),
    
    legend.title = element_text(
      face = "bold",
      size = 13
    ),
    
    legend.text = element_text(
      size = 9
    ),
    
    legend.position = "right",
    
    panel.grid = element_blank(),
    
    strip.background = element_rect(
      fill = "grey90",
      color = NA
    )
  ) +
  
  labs(
    title = "Top 20 Relative Abundant Species",
    x = "Sample ID",
    y = "Relative Abundance (%)",
    fill = "Genus_species"
  )

print(species_plot)

ggsave(
  file.path(output_dir, "Top20_Relative_Abundance_Species.png"),
  species_plot,
  width = 20,
  height = 10,
  dpi = 600,
  bg = "white"
)

write.csv(
  top20_species,
  file.path(output_dir, "Top20_Species.csv"),
  row.names = FALSE
)

# Alpha diversity by iron group

alpha_counts <- otu_tax |>
  select(Taxon, all_of(sample_columns)) |>
  group_by(Taxon) |>
  summarise(
    across(all_of(sample_columns), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  column_to_rownames("Taxon") |>
  as.matrix()

alpha_counts <- t(alpha_counts)
alpha_counts <- alpha_counts[
  rownames(alpha_counts) %in% metadata$SampleID, , drop = FALSE
]

alpha_meta <- metadata |>
  filter(
    SampleID %in% rownames(alpha_counts),
    !is.na(Iron)
  )

alpha_counts <- alpha_counts[alpha_meta$SampleID, , drop = FALSE]

iron_median <- median(alpha_meta$Iron, na.rm = TRUE)

alpha_meta <- alpha_meta |>
  mutate(
    Iron_Group = if_else(Iron <= iron_median, "Low Iron", "High Iron"),
    Iron_Group = factor(
      Iron_Group,
      levels = c("Low Iron", "High Iron")
    )
  )

alpha_values <- tibble(
  SampleID = rownames(alpha_counts),
  Shannon = vegan::diversity(alpha_counts, index = "shannon"),
  Simpson = vegan::diversity(alpha_counts, index = "simpson"),
  Observed = vegan::specnumber(alpha_counts)
) |>
  mutate(
    Pielou = if_else(
      Observed > 1,
      Shannon / log(Observed),
      NA_real_
    )
  ) |>
  left_join(
    alpha_meta |>
      select(SampleID, Iron, Iron_Group),
    by = "SampleID"
  )

alpha_long <- alpha_values |>
  pivot_longer(
    cols = c(Shannon, Simpson, Observed, Pielou),
    names_to = "Index",
    values_to = "Value"
  ) |>
  mutate(
    Index = factor(
      Index,
      levels = c("Shannon", "Simpson", "Observed", "Pielou"),
      labels = c(
        "Shannon diversity",
        "Simpson diversity",
        "Observed species",
        "Pielou evenness"
      )
    )
  )

alpha_wilcox <- alpha_long |>
  group_by(Index) |>
  group_modify(~ {
    tmp <- .x |> filter(!is.na(Value), !is.na(Iron_Group))
    n_low <- sum(tmp$Iron_Group == "Low Iron")
    n_high <- sum(tmp$Iron_Group == "High Iron")
    p_value <- if (n_low >= 1 && n_high >= 1) {
      wilcox.test(Value ~ Iron_Group, data = tmp, exact = FALSE)$p.value
    } else {
      NA_real_
    }
    tibble(n_low = n_low, n_high = n_high, p_value = p_value)
  }) |>
  ungroup() |>
  mutate(
    p_label = if_else(
      is.na(p_value),
      "p = NA",
      paste0("p = ", format.pval(p_value, digits = 2, eps = 0.001))
    )
  )

alpha_label_positions <- alpha_long |>
  group_by(Index) |>
  summarise(
    y_position = max(Value, na.rm = TRUE) +
      0.12 * (max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE) + 1e-8),
    .groups = "drop"
  ) |>
  left_join(alpha_wilcox, by = "Index")

alpha_plot_2A <- ggplot(
  alpha_long,
  aes(x = Iron_Group, y = Value, fill = Iron_Group)
) +
  geom_boxplot(
    width = 0.62,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  geom_jitter(
    width = 0.12,
    size = 1.8,
    alpha = 0.75
  ) +
  geom_text(
    data = alpha_label_positions,
    aes(x = 1.5, y = y_position, label = p_label),
    inherit.aes = FALSE,
    size = 3.4,
    fontface = "bold"
  ) +
  facet_wrap(~ Index, scales = "free_y", ncol = 2) +
  scale_fill_manual(
    values = c("Low Iron" = "#D95F5F", "High Iron" = "#4C9BCB")
  ) +
  labs(
    title = "Alpha Diversity by Iron Group",
    subtitle = paste0("Low/high iron split at median = ", signif(iron_median, 4)),
    x = NULL,
    y = "Alpha-diversity index"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.position = "none"
  )

print(alpha_plot_2A)

ggsave(
  file.path(output_dir, "Alpha_diversity_iron_groups.tiff"),
  alpha_plot_2A,
  width = 8, height = 6.5, dpi = 600, compression = "lzw"
)

write.csv(
  alpha_wilcox,
  file.path(output_dir, "Alpha_diversity_Wilcoxon_results.csv"),
  row.names = FALSE
)

# Spearman correlation with iron

iron_taxon_long <- otu_tax |>
  select(Taxon, all_of(sample_columns)) |>
  pivot_longer(
    cols = all_of(sample_columns),
    names_to = "SampleID",
    values_to = "Count"
  ) |>
  group_by(SampleID) |>
  mutate(
    Total_reads = sum(Count, na.rm = TRUE),
    Relative_Abundance = if_else(
      Total_reads > 0,
      Count / Total_reads,
      0
    )
  ) |>
  ungroup() |>
  group_by(Taxon, SampleID) |>
  summarise(
    Relative_Abundance = sum(Relative_Abundance, na.rm = TRUE),
    .groups = "drop"
  )

iron_metadata <- metadata |>
  filter(
    SampleID %in% sample_columns,
    !is.na(Iron),
    Timepoint %in% c("0 hours", "6 hours", "12 hours")
  ) |>
  mutate(
    Timepoint = factor(
      Timepoint,
      levels = c("0 hours", "6 hours", "12 hours")
    )
  )

top20_iron <- iron_taxon_long |>
  inner_join(
    iron_metadata |>
      select(SampleID),
    by = "SampleID"
  ) |>
  group_by(Taxon) |>
  summarise(
    Mean_Relative_Abundance = mean(Relative_Abundance, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(Mean_Relative_Abundance)) |>
  slice_head(n = 20)

iron_cor_data <- iron_taxon_long |>
  filter(Taxon %in% top20_iron$Taxon) |>
  inner_join(
    iron_metadata |>
      select(SampleID, Iron, Timepoint),
    by = "SampleID"
  )

spearman_iron <- iron_cor_data |>
  group_by(Taxon, Timepoint) |>
  group_modify(~ {
    dat <- .x |>
      filter(
        is.finite(Iron),
        is.finite(Relative_Abundance)
      )

    if (
      nrow(dat) < 3 ||
      dplyr::n_distinct(dat$Iron) < 2 ||
      dplyr::n_distinct(dat$Relative_Abundance) < 2
    ) {
      return(tibble(
        n = nrow(dat),
        rho = NA_real_,
        p_value = NA_real_
      ))
    }

    test <- suppressWarnings(
      cor.test(
        dat$Relative_Abundance,
        dat$Iron,
        method = "spearman",
        exact = FALSE
      )
    )

    tibble(
      n = nrow(dat),
      rho = unname(test$estimate),
      p_value = test$p.value
    )
  }) |>
  ungroup() |>
  mutate(
    significance = case_when(
      is.na(p_value) ~ "",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

taxon_order <- top20_iron |>
  arrange(Mean_Relative_Abundance) |>
  pull(Taxon)

spearman_iron <- spearman_iron |>
  mutate(
    Taxon = factor(Taxon, levels = taxon_order),
    Timepoint = factor(
      Timepoint,
      levels = c("0 hours", "6 hours", "12 hours")
    )
  )

iron_heatmap_2B <- ggplot(
  spearman_iron,
  aes(x = Timepoint, y = Taxon, fill = rho)
) +
  geom_tile(
    color = "black",
    linewidth = 0.45
  ) +
  geom_text(
    aes(label = significance),
    size = 4.5,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_gradient2(
    low = "#F44336",
    mid = "white",
    high = "#3B20E8",
    midpoint = 0,
    limits = c(-1, 1),
    breaks = c(-1, -0.5, 0, 0.5, 1),
    na.value = "grey90",
    name = expression("Spearman " * rho)
  ) +
  labs(
    title = "Spearman Correlation: Top 20 Species vs Iron",
    x = "Timepoint",
    y = "Bacterial species"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(
      face = "italic",
      color = "black",
      size = 9
    ),
    panel.grid = element_blank(),
    legend.title = element_text(face = "bold"),
    legend.key.height = grid::unit(1.2, "cm")
  )

print(iron_heatmap_2B)

ggsave(
  file.path(output_dir, "Top20_Spearman_Iron_heatmap.tiff"),
  iron_heatmap_2B,
  width = 7, height = 8, dpi = 600, compression = "lzw"
)

write.csv(
  spearman_iron,
  file.path(output_dir, "Spearman_Iron_results.csv"),
  row.names = FALSE
)

# Network topology

file_path <- network_file
network_data <- read_excel(file_path, sheet = "Sheet2")

print(network_data)
print(colnames(network_data))

network_data <- network_data %>%
  rename(
    Avg_path_length = `Avg.path length`, Avg_Degree = `Avg.Degree`)

network_data <- network_data %>%
  mutate(
    Time_clean = case_when(
      
      grepl(
        "^0|0\\s*hours?|0\\s*hrs?",
        as.character(Timepoints),
        ignore.case = TRUE) ~ "0 hours",
      
      grepl(
        "^6|6\\s*hours?|6\\s*hrs?",
        as.character(Timepoints),
        ignore.case = TRUE) ~ "6 hours",
      
      grepl(
        "^12|12\\s*hours?|12\\s*hrs?",
        as.character(Timepoints),
        ignore.case = TRUE) ~ "12 hours",
      
      TRUE ~ NA_character_))

print(unique(network_data$Timepoints))
print(unique(network_data$Time_clean))

network_data <- network_data %>%
  filter(
    Time_clean %in% c(
      "0 hours",
      "6 hours",
      "12 hours"))

network_data$Time_clean <- factor(
  network_data$Time_clean,
  levels = c(
    "0 hours",
    "6 hours",
    "12 hours"))

network_data <- network_data %>%
  mutate(
    Nodes = as.numeric(Nodes),
    Edges = as.numeric(Edges),
    Community = as.numeric(Community),
    Modularity = as.numeric(Modularity),
    Avg_path_length = as.numeric(Avg_path_length),
    Avg_Degree = as.numeric(Avg_Degree))

print(network_data)

print(
  colSums(
    is.na(
      network_data[
        c(
          "Nodes",
          "Edges",
          "Community",
          "Modularity",
          "Avg_path_length",
          "Avg_Degree")])))

network_long <- network_data %>%
  pivot_longer(
    
    cols = c(
      Nodes,
      Edges,
      Community,
      Modularity,
      Avg_path_length,
      Avg_Degree),
    names_to = "Metric",
    values_to = "Value")

network_long$Metric <- factor(
  network_long$Metric,
  levels = c(
    "Nodes",
    "Edges",
    "Community",
    "Modularity",
    "Avg_path_length",
    "Avg_Degree"),
  
  labels = c(
    "Nodes",
    "Edges",
    "Communities",
    "Modularity",
    "Average Path Length",
    "Average Degree"))

network_plot <- ggplot(
  network_long,
  aes(x = Time_clean, y = Value, group = 1)) +
  
  geom_line(linewidth = 1.5, alpha = 0.85) +
  
  geom_point(
    size = 5,
    shape = 21,
    fill = "white",
    stroke = 1.3) +
  
  facet_wrap(
    ~ Metric,
    scales = "free_y",
    ncol = 3) +
  
  scale_x_discrete(
    drop = FALSE,
    limits = c(
      "0 hours",
      "6 hours",
      "12 hours") ) +
  
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_markdown(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.spacing = unit(1.2, "lines")) +
  
  labs(
    title = "Network Topology Dynamics",
    subtitle = "Panchamritam fermentation: 0, 6 and 12 hours",
    x = NULL,
    y = "Metric Value")

print(network_plot)

output_file <- file.path(output_dir, "Network_Topology.png")

ggsave(
  filename = output_file,
  plot = network_plot,
  width = 12,
  height = 8,
  dpi = 600,
  bg = "white"
)

# Xanthomonas correlation heatmap

meta_data <- read_excel(
  meta_file
)

otu_table_raw <- read_excel(
  otu_file,
  sheet = "filtered_otu"
)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

otu_df <- as.data.frame(otu_table_raw)

tax_cols <- c(
  "OTU_ID", "Kingdom", "Phylum", "Class",
  "Order", "Family", "Genus", "Species"
)

missing_tax <- setdiff(tax_cols, names(otu_df))

if (length(missing_tax) > 0) {
  stop(
    "Missing taxonomy columns: ",
    paste(missing_tax, collapse = ", "))}

sample_cols <- setdiff(names(otu_df), tax_cols)

otu_df[sample_cols] <- lapply(
  otu_df[sample_cols],
  function(x) as.numeric(as.character(x)))

meta_df <- as.data.frame(meta_data)

required_meta <- c("SampleID", "Subject", "Timepoint")

if (!all(required_meta %in% names(meta_df))) {
  stop("Metadata must contain SampleID, Subject, and Timepoint.")
}

meta_df$SampleID <- trimws(as.character(meta_df$SampleID))

meta_df$Timepoint <- str_trim(as.character(meta_df$Timepoint))

meta_df$Timepoint <- case_when(
  str_detect(meta_df$Timepoint, "^0\\s*hours?$") ~ "0 hours",
  str_detect(meta_df$Timepoint, "^6\\s*hours?$") ~ "6 hours",
  str_detect(meta_df$Timepoint, "^12\\s*hours?$") ~ "12 hours",
  TRUE ~ meta_df$Timepoint)

timepoint_order <- c("0 hours", "6 hours", "12 hours")

meta_df <- meta_df %>%
  filter(Timepoint %in% timepoint_order) %>%
  distinct(SampleID, .keep_all = TRUE)

common_samples <- intersect(sample_cols, meta_df$SampleID)

if (length(common_samples) == 0) {
  stop("No matching sample IDs found between OTU table and metadata.")
}

message("Number of matched samples: ", length(common_samples))

meta_df <- meta_df %>%
  filter(SampleID %in% common_samples) %>%
  mutate(Timepoint = factor(Timepoint, levels = timepoint_order))

otu_df <- otu_df %>%
  mutate(
    Genus = str_trim(as.character(Genus)),
    Species = str_trim(as.character(Species)),
    
    Genus = na_if(Genus, ""),
    Genus = na_if(Genus, "NA"),
    Genus = na_if(Genus, "Unclassified"),
    Genus = na_if(Genus, "unclassified"),
    
    Species = if_else(
      is.na(Species) |
        Species == "" |
        Species %in% c("NA", "Unclassified", "unclassified"),
      "unclassified", Species),
    
    Taxon = if_else(
      is.na(Genus),
      NA_character_,
      paste(Genus, Species, sep = "_")
    )) %>%
  filter(!is.na(Taxon))

head(otu_df[, c("Genus", "Species", "Taxon")])

target_candidates <- unique(
  otu_df$Taxon[
    str_detect(
      otu_df$Taxon,
      regex("^Xanthomonas_", ignore_case = TRUE)
    )
  ]
)

cat("\nXanthomonas taxa found:\n")
print(target_candidates)

taxon_counts <- otu_df %>%
  select(Taxon, all_of(common_samples)) %>%
  group_by(Taxon) %>%
  summarise(
    across(all_of(common_samples), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )

count_matrix <- as.matrix(
  taxon_counts[, common_samples, drop = FALSE]
)

rownames(count_matrix) <- taxon_counts$Taxon

storage.mode(count_matrix) <- "numeric"

if (any(count_matrix < 0, na.rm = TRUE)) {
  stop("Negative OTU counts detected. Check the input data.")
}

count_matrix[is.na(count_matrix)] <- 0

sample_totals <- colSums(count_matrix)

if (any(sample_totals == 0)) {
  warning("Zero-total samples detected; they will be excluded.")
}

valid_samples <- names(sample_totals)[sample_totals > 0]

count_matrix <- count_matrix[, valid_samples, drop = FALSE]

relative_abundance <- sweep(
  count_matrix,
  MARGIN = 2,
  STATS = colSums(count_matrix),
  FUN = "/"
)

all_results <- list()

for (tp in timepoint_order) {
  
  cat("\nProcessing timepoint:", tp, "\n")
  
  tp_samples <- meta_df %>%
    filter(
      Timepoint == tp,
      SampleID %in% valid_samples
    ) %>%
    pull(SampleID)
  
  if (length(tp_samples) < 4) {
    warning("Too few samples at ", tp, "; skipping this timepoint.")
    next
  }
  
  if (!target_taxon %in% rownames(relative_abundance)) {
    stop("Target taxon missing from the abundance matrix.")
  }
  
  target_values <- relative_abundance[
    target_taxon, tp_samples, drop = TRUE
  ]
  
  if (length(unique(target_values)) < 2) {
    warning(
      "Target has no abundance variation at ", tp,
      "; correlations cannot be calculated."
    )
    next
  }
  
  other_taxa <- setdiff(
    rownames(relative_abundance),
    target_taxon
  )
  
  prevalence <- rowMeans(
    relative_abundance[
      other_taxa, tp_samples, drop = FALSE
    ] > 0
  )
  
  eligible_taxa <- other_taxa[prevalence > 0.30]
  
  cat("Samples:", length(tp_samples), "\n")
  cat("Taxa passing prevalence > 0.30:", length(eligible_taxa), "\n")
  
  if (length(eligible_taxa) == 0) next
  
  tp_results <- lapply(eligible_taxa, function(taxon) {
    
    x <- as.numeric(
      relative_abundance[taxon, tp_samples]
    )
    
    y <- as.numeric(target_values)
    
    keep <- is.finite(x) & is.finite(y)
    
    x <- x[keep]
    y <- y[keep]
    
    if (length(x) < 4 ||
        length(unique(x)) < 2 ||
        length(unique(y)) < 2) {
      return(NULL)
    }
    
    test <- suppressWarnings(
      cor.test(
        x,
        y,
        method = "spearman",
        exact = FALSE))
    
    data.frame(
      Timepoint = tp,
      Target = target_taxon,
      Bacterium = taxon,
      N = length(x),
      Prevalence = mean(x > 0),
      Spearman_rho = unname(test$estimate),
      P_value = test$p.value,
      stringsAsFactors = FALSE
    )
  })
  
  tp_results <- bind_rows(tp_results)
  
  if (nrow(tp_results) > 0) {
    all_results[[tp]] <- tp_results
  }
}

if (length(all_results) == 0) {
  stop(
    "No correlations could be calculated. Check target abundance, ",
    "sample numbers, and taxonomic labels."
  )
}

correlation_results <- bind_rows(all_results) %>%
  mutate(
    Timepoint = factor(Timepoint, levels = timepoint_order),
    Direction = if_else(
      Spearman_rho > 0,
      "Positive",
      "Negative"
    )
  ) %>%
  arrange(Timepoint, P_value)

write_csv(
  correlation_results,
  file.path(output_dir, "All_Spearman_Correlations.csv")
)

significant_results <- correlation_results %>%
  filter(P_value < 0.05) %>%
  arrange(Timepoint, P_value)

write_csv(
  significant_results,
  file.path(output_dir, "Significant_Spearman_p_below_0.05.csv")
)

print(significant_results)

if (nrow(significant_results) > 0) {
  
  heatmap_data <- significant_results %>%
    mutate(
      Bacterium = as.character(Bacterium),
      Timepoint = factor(
        as.character(Timepoint),
        levels = timepoint_order
      )
    )
  
  taxon_order <- heatmap_data %>%
    group_by(Bacterium) %>%
    summarise(
      max_abs_rho = max(abs(Spearman_rho)),
      .groups = "drop"
    ) %>%
    arrange(max_abs_rho) %>%
    pull(Bacterium)
  
  heatmap_data <- heatmap_data %>%
    mutate(
      Bacterium = factor(Bacterium, levels = taxon_order)
    )
  
  p_heatmap <- ggplot(
    heatmap_data,
    aes(
      x = Timepoint,
      y = Bacterium,
      fill = Spearman_rho
    )
  ) +
    geom_tile(
      color = "white",
      linewidth = 0.4
    ) +
    scale_fill_gradient2(
      low = "#D73027",      # Negative = red
      mid = "white",
      high = "#2166AC",     # Positive = blue
      midpoint = 0,
      limits = c(-1, 1),
      name = "Spearman\nrho"
    ) +
    labs(
      title = paste(
        "Spearman Correlation:",
        target_taxon,
        "vs Other Bacteria"
      ),
      subtitle = " p < 0.05; prevalence > 30%",
      x = "Timepoint",
      y = "Bacteria"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5
      ),
      axis.text.y = element_text(
        size = 8,
        color = "black"
      ),
      axis.title = element_text(face = "bold"),
      panel.grid = element_blank(),
      legend.title = element_text(face = "bold")
    )
  
  print(p_heatmap)
  
  ggsave(
    filename = file.path(
      output_dir,
      "Xanthomonas_Spearman_Heatmap.png"
    ),
    plot = p_heatmap,
    width = 9,
    height = max(5, 0.28 * length(taxon_order) + 2),
    units = "in",
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(output_dir, "Xanthomonas_Spearman_Heatmap.pdf"),
    plot = p_heatmap,
    width = 9,
    height = max(5, 0.28 * length(taxon_order) + 2),
    limitsize = FALSE)
  } else {
  message( "No correlations had unadjusted p < 0.05. ",
    "The heatmap was not generated.")}

for (tp in timepoint_order) {
  
  tp_significant <- significant_results %>%
    filter(as.character(Timepoint) == tp)
  
  safe_tp <- str_replace_all(tp, " ", "_")
  
  write_csv(
    tp_significant,
    file.path(output_dir, paste0("Significant_Correlations_", safe_tp, ".csv")))}

