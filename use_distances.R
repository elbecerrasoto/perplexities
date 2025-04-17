library(tidyverse)
library(dbscan)
library(Rtsne)
library(irlba) # dependendy of Rtsne
library(ggthemes)
library(glue)

system("trash *.pdf", intern = TRUE)

# Globals ----

RECALCULATE <- FALSE

MIN_PTS <- 4
PERPLEXITIES <- 2^c(3:8)
THREADS <- 12

SEED <- 424242
set.seed(SEED)

HITS_CHAR <- "ARCH.tsv"

# Helpers ----

plot_tsne <- function(plot_data, ititle, isubtitle) {
  mysam <- function(x) {
    if_else(
      sample(c(TRUE, FALSE), length(x), replace = TRUE, prob = c(0.32, 0.68)),
      x,
      NA
    )
  }

  ggplot(plot_data) +
    geom_point(aes(x = mysam(V1), y = mysam(V2), color = fct_lump(cluster, 8)), size = 5, shape = 1) +
    geom_point(aes(x = V1, y = V2), alpha = 1 / 8, size = 0.64) +
    theme_fivethirtyeight(base_size = 18) +
    guides(color = guide_legend(title = NULL)) +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    ) +
    scale_color_manual(values = paletteer::paletteer_d("ggthemes::few_Light")) +
    labs(
      title = ititle,
      subtitle = isubtitle,
      caption = "author: Becerra-Soto E.",
      x = "",
      y = ""
    )
}

# Main ----

hits_char <- read_tsv(HITS_CHAR)

wdl <- read_rds("wdl.Rds")
HDB_wdl <- hdbscan(as.dist(wdl), minPts = MIN_PTS)

x <- as_tibble(matrix(wdl, nrow = nrow(wdl), ncol = ncol(wdl)))
x$noise <- rnorm(nrow(wdl), mean = 0.0, sd = 1e-5)

perplexity_search <- vector(mode = "list", length = length(PERPLEXITIES))
perplexity_names <- str_c("perp_", PERPLEXITIES)
names(perplexity_search) <- perplexity_names

if (RECALCULATE) {
  for (i in seq_along(PERPLEXITIES)) {
    tsne <- Rtsne(x,
      dims = 2, perplexity = PERPLEXITIES[i], verbose = TRUE,
      partial_pca = TRUE, num_threads = THREADS
    )
    perplexity_search[[i]] <- tsne
  }

  write_rds(perplexity_search, "perplexity_search.Rds")
} else {
  perplexity_search <- read_rds("perplexity_search.Rds")
}


# 1111111111111111111

perplexity_plots <- vector(mode = "list", length = length(PERPLEXITIES))
names(perplexity_plots) <- perplexity_names

for (i in seq_along(perplexity_search)) {
  tsne <- perplexity_search[[i]]
  plot_data <- as_tibble(tsne$Y)
  plot_data$neID <- rownames(wdl)
  plot_data$cluster <- as_factor(HDB_wdl$cluster)
  plot_data$outlier <- HDB_wdl$cluster == 0

  title <- glue("HDB {names(perplexity_search)[i]}")
  subtitle <- "HDB_clusters"
  perplexity_plots[[i]] <- plot_tsne(plot_data, title, subtitle)
}


# 2222222222222222

curated_plots <- vector(mode = "list", length = length(PERPLEXITIES))
names(curated_plots) <- perplexity_names

for (i in seq_along(perplexity_search)) {
  tsne <- perplexity_search[[i]]
  plot_data <- as_tibble(tsne$Y)
  plot_data$neID <- rownames(wdl)

  plot_data <- left_join(plot_data, hits_char, join_by(neID))
  plot_data$cluster <- as_factor(plot_data$curated)

  title <- glue("PFAM {names(curated_plots)[i]}")
  subtitle <- "PFAM curated clusters"
  curated_plots[[i]] <- plot_tsne(plot_data, title, subtitle)
}

# 33333333333333333333333333

phylum_plots <- vector(mode = "list", length = length(PERPLEXITIES))
names(phylum_plots) <- perplexity_names

for (i in seq_along(perplexity_search)) {
  tsne <- perplexity_search[[i]]
  plot_data <- as_tibble(tsne$Y)
  plot_data$neID <- rownames(wdl)

  plot_data <- left_join(plot_data, hits_char, join_by(neID))
  plot_data$cluster <- as_factor(plot_data$phylum)

  title <- glue("Phylum {names(phylum_plots)[i]}")
  subtitle <- "Phylum clusters"
  phylum_plots[[i]] <- plot_tsne(plot_data, title, subtitle)
}

# Plot

for (name_perp_plot in names(perplexity_plots)) {
  ggsave(glue("hdbclust_{name_perp_plot}.pdf"), plot = perplexity_plots[[name_perp_plot]], width = 8.5, height = 11, units = "in", dpi = 300)
}

for (name_phyl_plot in names(phylum_plots)) {
  ggsave(glue("phylum_{name_phyl_plot}.pdf"), plot = phylum_plots[[name_perp_plot]], width = 8.5, height = 11, units = "in", dpi = 300)
}

for (name_curat_plot in names(curated_plots)) {
  ggsave(glue("curated_{name_curat_plot}.pdf"), plot = curated_plots[[name_perp_plot]], width = 8.5, height = 11, units = "in", dpi = 300)
}

# combine
system("pdfunite hdbclust_perp_* phylum_perp_* curated_perp_* perplexities.pdf", intern = TRUE)
