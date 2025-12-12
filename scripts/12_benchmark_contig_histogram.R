library(ggplot2)
library(Biostrings)
library(purrr)
library(tibble)
library(dplyr)

#----------------SUMMARY OF ALL ISOLATES------------
data_dir <- "/g/typas/Personal_Folders/Neeka/Model/data/08_benchmark_simulator/ERP105624_selection/assembled_genomes"

#need to load in contig file from each assembly (.fasta)
#then detect all lines that are node (aka start of one contig)
#and then count the number of characters until next node (aka length of one contig)

summary_files <- Sys.glob(file.path(data_dir, "*/contigs.fasta"))

length(summary_files)
summary_files

node_df <- map_dfr(summary_files, function(f) {
  s <- readDNAStringSet(f)
  tibble(
    sample = basename(dirname(f)),
    node   = names(s),
    length = width(s)
  )
})

p <- ggplot(node_df, aes(x=length)
  ) + 
  geom_histogram(
    bins = 50,                  # adjust number of bins (more = finer detail)
    fill = "#2E86AB",           # soft blue fill
    color = "white"             # white outline between bins
  ) +
   scale_x_log10(
   ) +
  labs(
    title = "Contig length distribution",
    x = "Contig length (bp)",
    y = "Number of contigs"
  ) +
  theme_minimal(
  ) 

ggsave(
  filename = file.path(data_dir, "contigs_histogram.png"),
  plot = p,
  width = 10,
  height = 7,
  dpi = 300
)

node_df %>% 
  group_by(sample) %>% 
  summarise(n_contigs = n(),
            total_bp = sum(length),
            N50 = quantile(length, 0.5))
#the contig length such that 50% of the total assembly length (in base pairs) is contained in contigs of this length or longer
#a measure of assembly contiguity: larger values mean fewer, longer contigs.

#----------------ISOLATES ALONE------------
# Get all immediate subdirectories (one per species)
species_dirs <- list.dirs(data_dir, full.names = TRUE, recursive = FALSE)

# Initialize an empty data frame
contig_lengths <- tibble()

for (sp_dir in species_dirs) {
  contig_file <- file.path(sp_dir, "contigs.fasta")
  
  if (file.exists(contig_file)) {
    s <- readDNAStringSet(contig_file)
    contig_lengths <- bind_rows(contig_lengths, tibble(
      species = basename(sp_dir),
      length = width(s)
    ))
  }
}

# 📊 Plot a histogram of contig lengths, faceted per species
p <- ggplot(contig_lengths, aes(x = length)) + 
  geom_histogram(bins = 50, fill = "#2E86AB", color = "white") +
  scale_x_log10() +
  labs(
    title = "Contig length distribution per species",
    x = "Contig length (bp, log scale)",
    y = "Number of contigs"
  ) +
  theme_minimal() +
  facet_wrap(~species, nrow = 2, ncol = 5)

ggsave(
  filename = file.path(data_dir, "contigs_histogram_by_species.png"),
  plot = p,
  width = 12,
  height = 8,
  dpi = 300
) 