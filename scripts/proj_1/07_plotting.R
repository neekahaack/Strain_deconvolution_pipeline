library(ggplot2)
library(dplyr)
library(readr)
library(forcats)
library(tidyr)
library(stringr)
library(purrr)

#need to pivot long like in 06 script, but do it here instead

taxonomy_base <- "/g/typas/Personal_Folders/Neeka/Model/data/02_taxonomy/"
matrices_dir  <- "/g/typas/Personal_Folders/Neeka/Model/data/05_aligned/pairwise_by_genus"
pivoted_matrix <- "/g/typas/Personal_Folders/Neeka/Model/data/06_analysis/pairwise_identity_by_genus_16s_pivoted_long.csv"

# -------- 1) Load taxonomy and extract genus --------
summary_files <- list.files(
  taxonomy_base,
  pattern = "gtdbtk\\.bac120\\.summary\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

taxonomy <- map_dfr(summary_files, ~ read_tsv(.x, show_col_types = FALSE))
# expect at least columns: user_genome, classification
taxonomy <- taxonomy %>%
  mutate(genus = str_match(classification, "g__([^;]+)")[,2]) %>%
  select(user_genome, genus)

# -------- 2) Read every per-genus identity matrix and pivot to long --------
matrix_files <- list.files(
  matrices_dir,
  pattern = "\\.pairwise_identity\\.csv$",
  full.names = TRUE
)

if (length(matrix_files) == 0) {
  stop("No per-genus matrix CSVs found in matrices_dir")
}

read_one_matrix_long <- function(path) {
  df <- read_csv(path, show_col_types = FALSE)
  # first column is the row IDs (sequence IDs). Make sure it's named 'seq_i'
  first_col <- names(df)[1]
  df <- df %>% rename(seq_i = !!first_col)

  # pivot wide -> long
  df_long <- df %>%
    pivot_longer(
      cols = -seq_i,
      names_to = "seq_j",
      values_to = "identity"
    )

  # add genus label from filename (e.g., "Bacteroides.pairwise_identity.csv" -> "Bacteroides")
  genus_from_file <- basename(path) %>% sub("\\.pairwise_identity\\.csv$", "", .)
  df_long$genus_matrix <- genus_from_file
  df_long
}

long_df <- map_dfr(matrix_files, read_one_matrix_long)

# -------- 3) Clean and annotate --------
# identity as numeric; drop NAs and self-comparisons
long_df <- long_df %>%
  mutate(identity = suppressWarnings(as.numeric(identity))) %>%
  filter(!is.na(identity), seq_i != seq_j)

# extract user_genome from sequence IDs (token after last "_____")
long_df <- long_df %>%
  mutate(
    user_genome_i = sub(".*_____", "", seq_i),
    user_genome_j = sub(".*_____", "", seq_j)
  ) %>%

  # add genus_i and genus_j from taxonomy
  left_join(taxonomy %>% rename(user_genome_i = user_genome, genus_i = genus),
            by = "user_genome_i") %>%
  left_join(taxonomy %>% rename(user_genome_j = user_genome, genus_j = genus),
            by = "user_genome_j")

# keep upper triangle only (stable ordering by string)
long_df <- long_df %>%
  filter(seq_i < seq_j)%>%
# drop rows with missing genus (should be rare)
  filter(!is.na(genus_i), !is.na(genus_j))%>%
# identities should be 0..100; keep plausible range
  filter(identity >= 0, identity <= 100)

# -------- 4) Save combined long table --------
write_csv(long_df, pivoted_matrix)

# Optional: quick peek
print(glue::glue("Rows: {nrow(long_df)}   Genera (from files): {dplyr::n_distinct(long_df$genus_matrix)}"))
print(head(long_df, 5))

plot_0 <- ggplot(
  long_df,
  aes(
    x = fct_reorder(genus_i, identity, .fun = median),#reorder genus_i by median value of identity
    y = identity,
    fill = genus_i)
  ) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8
  ) +
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Genus", y = "Pairwise Identity (%)")

ggsave(
  filename = "/g/typas/Personal_Folders/Neeka/Model/data/07_plotting/by_genus_boxplot_redo.png",
  plot = plot_0,
  width = 10,
  height = 14,
  dpi = 300
)
  
#----------------intragenus boxplot----------------

data <- read_csv("/g/typas/Personal_Folders/Neeka/Model/data/06_analysis/pairwise_identity_16s_pivoted_long.csv")
#but also indirectly intergenus, by comparing their within-genus variability and medians.
intragenus <- data %>%
  filter(genus_i == genus_j) 

# Count number of comparisons per genus
counts <- intragenus %>%
  group_by(genus_i) %>%
  summarize(n = n())

nrow(counts)    #distinct genera
sum(counts$n)   #intersections counted

plot <- ggplot(
  intragenus,
  aes(
    x = fct_reorder(genus_i, identity, .fun = median),#reorder genus_i by median value of identity
    y = identity,
    fill = genus_i)
  ) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8
  ) +
  geom_text(
    data = counts,
    aes(x = genus_i, y = 102, label = paste0("n=", n)),  # places counts above boxes
    size = 2.8,
    hjust = 0
  ) +
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Genus", y = "Pairwise Identity (%)")

ggsave(
  filename = "/g/typas/Personal_Folders/Neeka/Model/data/07_plotting/intragenus_boxplot.png",
  plot = plot,
  width = 10,
  height = 14,
  dpi = 300
)

#----------------intergenus boxplot----------------
intergenus <- data %>%
  filter(genus_i != genus_j) 

#count number of comparisons per genus
counts <- intergenus %>%
  group_by(genus_i) %>%
  summarize(n = n())

nrow(counts)    #is like 188 distinct genera
sum(counts$n)   #is number of combinations/intersections/identity values being counted

plot_2<- ggplot(
  intergenus,
  aes(
    x = fct_reorder(genus_i, identity, .fun = median),#reorder genus_i by median value of identity
    y = identity,
    fill = genus_i)
  ) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8
  ) +
  geom_text(
    data = counts,
    aes(x = genus_i, y = 102, label = paste0("n=", n)),  # places counts above boxes
    size = 2.8,
    hjust = 0
  ) +
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Genus", y = "Pairwise Identity (%)")

ggsave(
  filename = "/g/typas/Personal_Folders/Neeka/Model/data/07_plotting/intergenus_boxplot.png",
  plot = plot_2,
  width = 10,
  height = 18,
  dpi = 300
)

#--------------intergenus heatmap----------------

intergenus_summary <- intergenus %>%
  group_by(genus_i, genus_j) %>%
  summarize(median_identity = median(identity), .groups = "drop")

plot_3 <- ggplot(intergenus_summary, aes(x = genus_i, y = genus_j, fill = median_identity)) +
  geom_tile() +
  coord_equal() +
  scale_fill_viridis_c() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))  #rotate x labels

ggsave(
  filename = "/g/typas/Personal_Folders/Neeka/Model/data/07_plotting/intergenus_heatmap.png",
  plot = plot_3,
  width = 25,
  height = 25,
  dpi = 300
)
