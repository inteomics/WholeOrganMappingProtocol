# Load data
set.seed(1918)

library(ComplexHeatmap)
library(patchwork)
library(phangorn)
library(ape)
library(e1071)
library(ade4)
library(tidygraph)
library(ggraph)
library(tidyverse)

samples_data <- readRDS("../data/samples_data.Rds") %>%
  filter(map == "map19")
ur_samples <- filter(samples_data, group != "MUSCLE")
filtered_silent_nonsilent <- readRDS("../data/filtered_mutations_data.map19.Rds")
VAF_mat <- readRDS("../data/VAF_mat.map19.Rds")

dim(VAF_mat)

# Binary martix prepareation

VAF_mat <- VAF_mat[filtered_silent_nonsilent$mut_id, ur_samples$sam_id]
colnames(VAF_mat) <- str_replace(colnames(VAF_mat), "map19_", "")

calc_nonzero_fract <- function(t) {
  sum(VAF_mat > t) / (nrow(VAF_mat) * ncol(VAF_mat))
}

treshold_stats <- tibble(
  treshold = seq(0.0001, 0.1, length.out = 2000),
  fraction1 = map_dbl(treshold, calc_nonzero_fract)
)

ggplot(treshold_stats, aes(treshold, fraction1)) +
  geom_point() +
  scale_x_log10() +
  labs(y = "% of non-zero elements of mutations' matrix")

treshold <- 0.01
VAF_mat_bin <- matrix(0, nrow = nrow(VAF_mat), ncol = ncol(VAF_mat))
colnames(VAF_mat_bin) <- colnames(VAF_mat)
rownames(VAF_mat_bin) <- rownames(VAF_mat)
VAF_mat_bin[VAF_mat > treshold] <- 1

cat("Selected solution:")
sum(VAF_mat_bin) / (nrow(VAF_mat_bin) * ncol(VAF_mat_bin))

# Order of columns/rows in this matrix affects results. Use order from the first
# analysis with 3 branches, based on old names
col_order <- c("1-D8","12-E14", "13-F7",  "18-F13", "19-F14", "2-D10",  "20-G7",
  "27-G14", "28-H7",  "29-H8",  "3-D11", "32-H11", "33-H12", "34-H13", "35-H14",
  "36-I9",  "37-I12", "4-D12",  "5-E7", "6-E8", "7-E9", "10-E12", "11-E13",
  "17-F12", "24-G11", "25-G12", "26-G13", "30-H9", "31-H10", "9-E11", "14-F8",
  "15-F9", "16-F11", "21-G8", "22-G9", "23-G10", "8-E10") %>%
  str_match("(?<=\\-)[:graph:]+")
VAF_mat_bin <- VAF_mat_bin[, col_order]

# Trees construction

hamming_dist <- hamming.distance(t(VAF_mat_bin))

sample_classification_data <- samples_data %>%
    select(region, group) %>%
    deframe() %>%
    .[rownames(hamming_dist)]

right_anno <- rowAnnotation(
  Group = sample_classification_data,
  col = list(
    Group = c(`NU/LGIN` = "#3dbce7", HGIN = "#ffac49", UC = "black")
  )
)

ht <- Heatmap(
  hamming_dist,
  right_annotation = right_anno
)


set.seed(1918)
trees <- list(
  UPGMA = upgma(hamming_dist),
  NJ = NJ(hamming_dist)
)
map(trees, is.rooted)

mutdat <- phyDat(t(VAF_mat_bin), type = "USER", levels = c(0, 1))

parsimony_trees <- readRDS(str_c("results/parsimony_trees_t=", treshold, ".Rds"))

identical(ur_samples$region, parsimony_trees$UPGMA$tip.label)
identical(ur_samples$region, parsimony_trees$NJ$tip.label)

tip_colors <- ur_samples %>%
  mutate(color = case_when(
    group == "NU/LGIN" ~ "#3dbce7",
    group == "HGIN" ~ "#ffac49",
    group == "UC" ~ "black",
  )) %>%
  pull(color)

iwalk(parsimony_trees, 
  ~plot(.x, "unrooted",
    main = str_c("MP unrooted tree (from ", .y, ")"),
    cex = .5, font = 2, no.margin = F, lab4ut = "axial",
    tip.color = tip_colors
  )
)


# Diagnose trees

dist_vec <- as.vector(hamming_dist)
UPGMA_dist_vec <- as.vector(cophenetic(parsimony_trees$UPGMA))
NJ_dist_vec <- as.vector(cophenetic(parsimony_trees$NJ))

qplot(dist_vec, UPGMA_dist_vec) +
  geom_abline(slope = 1, intercept = 0) +
  labs(
    title = "Is UPGMA appropriate?",
    x = "original pairwise distances",
    y = "pairwise distances on the tree"
  )

qplot(dist_vec, NJ_dist_vec) +
  geom_abline(slope = 1, intercept = 0) +
  labs(
    title = "Is NJ appropriate?",
    x = "original pairwise distances",
    y = "pairwise distances on the tree"
  )
# dev.off()

# Bootstraping

# my_boots <- boot.phylo(parsimony_trees$NJ, t(VAF_mat_bin),
#   function(x) {
#     hamming_dist <- hamming.distance(x)
#     tree_NJ <- NJ(hamming_dist)
#     tree_NJ
#     mutdat <- phyDat(x, type = "USER", levels = c(0, 1))
#     optim.parsimony(tree_NJ, mutdat)
#   },
#   B = 10
# )

bootstrap_file <- str_c("results/bootstrap_t=", treshold, ".Rds")

if (file.exists(bootstrap_file)) {
  bootstrap <- readRDS(str_c("results/bootstrap_t=", treshold, ".Rds"))
} else {
  bootstrap <- list()
    
  bootstrap$NJ <- bootstrap.phyDat(mutdat,
    function(x) {
      y <- x %>%
        as.data.frame() %>%
        map_df(as.numeric) %>%
        as.matrix()
      hamming_dist <- hamming.distance(t(y))
      tree_NJ <- NJ(hamming_dist)
      optim.parsimony(tree_NJ, x)
    }
  )
  
  bootstrap$NJ_unopt <- bootstrap.phyDat(mutdat,
    function(x) {
      y <- x %>%
        as.data.frame() %>%
        map_df(as.numeric) %>%
        as.matrix()
      hamming_dist <- hamming.distance(t(y))
      NJ(hamming_dist)
    }
  )
  
  bootstrap$UPGMA <- bootstrap.phyDat(mutdat,
    function(x) {
      y <- x %>%
        as.data.frame() %>%
        map_df(as.numeric) %>%
        as.matrix()
      hamming_dist <- hamming.distance(t(y))
      tree_upgma <- upgma(hamming_dist)
      optim.parsimony(unroot(tree_upgma), x)
    }
  )
  
  saveRDS(bootstrap, file = bootstrap_file)
}

plotBS(parsimony_trees$NJ, bootstrap$NJ, tip.color = tip_colors, main = "NJ")

plotBS(parsimony_trees$UPGMA, bootstrap$UPGMA, tip.color = tip_colors, main = "UPGMA")


branches <- openxlsx::read.xlsx("../map26 - resources/map26.tree_branches.xlsx")
colors <- branches %>%
  mutate(
    tip = str_replace(sam_id, "map26_", ""),
    color = case_when(
      branch == "eps" ~ "#1d457f",
      branch == "gamma" ~ "#cc5c76",
      branch == "delta" ~ "#f9ad2a"
    )
  ) %>%
  select(tip, color) %>%
  deframe()

pdf("results/9 bootstrapping trees.pdf", width = 20, height = 20)
plotBS(parsimony_trees$NJ, bootstrap$NJ, tip.color = colors[parsimony_trees$NJ$tip.label], main = "NJ")

patchwork::wrap_elements(~plot(bootstrap$NJ[[1]], "unrooted", rot = -25, tip.color = colors[bootstrap$NJ[[1]]$tip.label])) + 
patchwork::wrap_elements(~plot(bootstrap$NJ[[2]], "unrooted", rot = 200, tip.color = colors[bootstrap$NJ[[2]]$tip.label])) + 
patchwork::wrap_elements(~plot(bootstrap$NJ[[3]], "unrooted", rot = -80, tip.color = colors[bootstrap$NJ[[3]]$tip.label])) + 
patchwork::wrap_elements(~plot(bootstrap$NJ[[4]], "unrooted", rot = -25, tip.color = colors[bootstrap$NJ[[4]]$tip.label])) + 
patchwork::wrap_elements(~plot(bootstrap$NJ[[5]], "unrooted", rot = -75, tip.color = colors[bootstrap$NJ[[5]]$tip.label])) + 
patchwork::wrap_elements(~plot(bootstrap$NJ[[6]], "unrooted", rot = -25, tip.color = colors[bootstrap$NJ[[6]]$tip.label])) + 
patchwork::wrap_elements(~plot(bootstrap$NJ[[7]], "unrooted", rot = -120, tip.color = colors[bootstrap$NJ[[7]]$tip.label])) + 
patchwork::wrap_elements(~plot(bootstrap$NJ[[8]], "unrooted", rot = -150, tip.color = colors[bootstrap$NJ[[8]]$tip.label])) + 
patchwork::wrap_elements(~plot(bootstrap$NJ[[9]], "unrooted", rot = -90, tip.color = colors[bootstrap$NJ[[9]]$tip.label]))
dev.off()

#Custom visualizations

tidy_trees <- map(parsimony_trees, as_tbl_graph)
tidy_trees$NJ

tidy_trees$NJ %>%
  activate(nodes) %>%
  mutate(sam_id = name) %>%
  left_join(samples_data) %>%
  ggraph(layout = "tree") +
  geom_node_point(aes(filter = is.na(group))) +
  geom_node_label(aes(color = group, label = name, filter = !is.na(group)),
                  fill = "white", check_overlap = TRUE) +
  geom_edge_diagonal() +
  scale_color_manual(values = c("#3dbce7", "#ffac49", "black"), breaks = c("NU/LGIN", "HGIN", "UC"))
