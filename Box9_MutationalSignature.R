foldername <- "Analysis"
if (!dir.exists(foldername)) dir.create(foldername)
set.seed(4)
colors <- list(
    variant_spread = c(
        `31+ samples` = "red4",
        `21 - 30 samples` = "firebrick3",
        `11 - 20 samples` = "darkorange1",
        `2 - 10 samples` = "goldenrod1",
        `private` = "dodgerblue4"
    ),
    blca_3_class_red = c(`NU/LGIN` = "#fdbb84", HGIN = "#ef6548", UC = "#b40000"),
    blca_4_class_red = c(NU = "#fee8c8", LGIN = "#fdbb84", HGIN = "#ef6548", UC = "#b30000"),
    blca_3_class_orig = c(`NU/LGIN` = "#3dbce7", HGIN = "#ffac49", UC = "black"),
    blca_codes = c(N = "#edf8fb", D1 = "#bfd3e6", D2 = "#9ebcda", D3 = "#8c96c6", D4 = "#8856a7", D5 = "#810f7c"),
    mutational_landscape = c("#3952a4", "#eff0f1", "#ee1f25")
)

# ======================================================================
# ====================================== ANALYSES OF MUTATION SIGNATURES
# ======================================================================
set.seed(1918)
library(ComplexHeatmap)
library(patchwork)
library(mutSignatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(scales)
library(tidyverse)
hg38 <- BSgenome.Hsapiens.UCSC.hg38
theme_set(theme_minimal())
treshold <- 0.01
#   Read in data
samples_data <- read_tsv("samples_data.tsv") |>
    filter(patient_id == "M19") |>
    mutate(
        sample_id = str_replace(sample_id, "M19_", ""),
        group = parse_factor(group3, levels = c("NU/LGIN", "HGIN", "UC"))
    )
ur_samples <- filter(samples_data, group != "MUSCLE")
filtered_mutations_data <- readRDS("filtered_variants.Rds") |>
    filter(patient_id == "M19")
VAF_mat <- readRDS("VAF_matrix.Rds")
colnames(VAF_mat) <- colnames(VAF_mat) |>
    str_replace("map19_", "")
#   Filter data for only SNV mutations
nucleotides <- c("A", "T", "C", "G")
filtered_mutations_data <- filtered_mutations_data |>
    filter(
        effect %in% c("nonsilent", "null"),
        # variant_type == "SNV",
        ref %in% nucleotides,
        alt %in% nucleotides
    )
#   Join VAF data and mutations data
VAF_mat <- VAF_mat |>
    as.data.frame() |>
    rownames_to_column("mut_id")
VAF_mat$mut_id <- VAF_mat$mut_id %>%
    {
        gsub(":.*?:", "-", .)
    } %>%
    {
        gsub("(.)>", "-\\1>", .)
    } %>%
    {
        gsub(">", "-", .)
    }
VAF_mat <- VAF_mat %>%
    filter(mut_id %in% filtered_mutations_data$mutation_id)
#   Attach nucleotide context for mutations
mutations_data <- attachContext(filtered_mutations_data,
    chr_colName = "chrom",
    start_colName = "pos",
    end_colName = "pos",
    nucl_contextN = 3,
    BSGenomeDb = hg38
)
mutations_data <- removeMismatchMut(mutations_data,
    refMut_colName = "ref",
    context_colName = "context",
    refMut_format = "N"
)
mutations_data <- attachMutType(mutations_data,
    ref_colName = "ref",
    var_colName = "alt",
    context_colName = "context"
)
mutations_data <- mutations_data |>
    mutate(
        sig_substitution = str_sub(mutType, 3, 5) |>
            parse_factor(levels = c("T>G", "T>C", "T>A", "C>T", "C>G", "C>A"))
    )
VAF_mat <- VAF_mat |>
    as.data.frame() |>
    pivot_longer(-mut_id, names_to = "sam_id", values_to = "VAF")
VAF_tbl <- VAF_mat |>
    left_join(mutations_data %>% select(mutation_id, mut_type = mutType, sig_substitution), by = c("mut_id" = "mutation_id")) |>
    left_join(samples_data %>% select(sample_id, group), by = c("sam_id" = "sample_id"))
#---Fig a) - SNVs by substitution and region**
#   mutation considered present if VAF > 0.01. Results strongly depend on this threshold
substitution_colors <- c(
    `C>A` = "deepskyblue1", `C>G` = "black", `C>T` = "#ee1c25",
    `T>A` = "Gray", `T>C` = "#59ba46", `T>G` = "#f9bfcc"
)
substitution_counts <- VAF_tbl %>%
    filter(VAF > treshold) %>%
    group_by(sig_substitution, group) %>%
    count()
p <- list()
p$a <- ggplot(substitution_counts, aes(group, n, fill = sig_substitution)) +
    geom_bar(stat = "identity", position = "fill", width = 0.7) +
    scale_y_continuous(labels = percent_format()) +
    scale_fill_manual(values = substitution_colors) +
    theme_minimal() +
    theme(
        legend.position = "top"
    ) +
    labs(
        y = "SNVs",
        fill = "Type"
    )
p$a
ggsave(pathify("a_map19-stackedBarplotForSNVs.pdf"))
ggsave(paste0(foldername, "/[3A TOP] Stacked Barplot For SNVs.pdf"))
#   Proportion test
n_snv_sum <- substitution_counts %>%
    pivot_wider(names_from = sig_substitution, values_from = n) %>%
    column_to_rownames("group")
grade_total <- apply(n_snv_sum, 1, sum)
cat("Grade total\n")
grade_total
prop_p <- apply(n_snv_sum, 2, function(x) prop.test(x, grade_total)$p.value)
cat("\nprop_p\n")
prop_p
adj_prop_p <- p.adjust(prop_p, method = "holm")
cat("\nround(adj_prop_p, 3)\n")
round(adj_prop_p, 3)
pdat_perc <- sweep(n_snv_sum, 1, grade_total, "/")
pdat_perc
data.frame(t(pdat_perc), prop_p, adj_prop_p) %>%
    rownames_to_column("substitution") %>%
    write_tsv(paste0(foldername, "/[3A TOP] StackedBarplotProportionTest.tsv"))

#   b - the barplots of base substitution with flanking sequences**
sig_substitution_counts <- VAF_tbl %>%
    filter(VAF > treshold) %>%
    group_by(group, mut_type) %>%
    count() %>%
    group_by(group) %>%
    mutate(
        mut_perc = n / sum(n),
        substitution = str_sub(mut_type, 3, 5)
    )
#   plot 3B, flanking bases
p$b <- ggplot(sig_substitution_counts, aes(mut_type, mut_perc)) +
    geom_bar(stat = "identity", width = 0.2, color = "black") +
    scale_y_continuous(labels = percent_format()) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 4)
    ) +
    labs(
        y = "SNVs",
        x = "Type"
    ) +
    facet_grid(group ~ substitution, space = "free", scales = "free_x")
p$b
ggsave(paste0(foldername, "/[3B TOP] baseSubstitutionFlankingBasesPerc-barplot.pdf"))
#   Wilcoxon rank sum test
#   Prepare object the same as Hui's
all_mut_types <- expand_grid(
    prev = nucleotides,
    mut = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
    after = nucleotides
) %>%
    transmute(mut_type = str_c(prev, "[", mut, "]", after)) %>%
    pull(mut_type)
snv_nb_sam_h <- VAF_tbl %>%
    filter(VAF > treshold) %>%
    mutate(
        sam_id = parse_factor(sam_id, levels = ur_samples$sample_id),
        mut_type = parse_factor(mut_type, levels = all_mut_types)
    ) %>%
    group_by(sam_id, mut_type) %>%
    count() %>%
    ungroup() %>%
    complete(sam_id, mut_type, fill = list(n = 0)) %>%
    pivot_wider(id_cols = "sam_id", names_from = "mut_type", values_from = n)
this_map_sam <- ur_samples %>%
    dplyr::rename(type = group) %>%
    mutate(type = forcats::fct_drop(type))
mat.baseSubstitutionFlankingBases <- column_to_rownames(snv_nb_sam_h, "sam_id")
mat.baseSubstitutionFlankingBases <- mat.baseSubstitutionFlankingBases[this_map_sam$sample_id, ]
mat.baseSubstitutionFlankingBasesPerc <- sweep(
    mat.baseSubstitutionFlankingBases, 1,
    apply(mat.baseSubstitutionFlankingBases, 1, sum), "/"
)
kruskal_pvalue <- apply(mat.baseSubstitutionFlankingBasesPerc, 2, function(x) {
    tdat <- data.frame(perc = x, type = this_map_sam$type)
    kruskal.test(perc ~ type, tdat)$p.value
})
sort(round(kruskal_pvalue, 3))
type_list <- levels(this_map_sam$type)
wilcox_pvalue <- matrix(NA, nrow = 3, ncol = ncol(mat.baseSubstitutionFlankingBasesPerc))
colnames(wilcox_pvalue) <- colnames(mat.baseSubstitutionFlankingBasesPerc)
rownames(wilcox_pvalue) <- 1:3
k <- 1
for (i in 1:(length(type_list) - 1)) {
    for (j in (i + 1):length(type_list)) {
        idx <- this_map_sam$type %in% c(type_list[i], type_list[j])
        pvalue <- apply(mat.baseSubstitutionFlankingBasesPerc[idx, ], 2, function(x) {
            tdat <- data.frame(perc = x, type = this_map_sam$type[idx])
            tdat$type <- factor(tdat$type)
            wilcox.test(perc ~ type, tdat)$p.value
        })

        wilcox_pvalue[k, ] <- pvalue
        rownames(wilcox_pvalue)[k] <- paste0(type_list[i], " vs ", type_list[j])
        k <- k + 1
    }
}
### correct multipling testing for all tests
wilcox_adj_pvalue <- matrix(p.adjust(wilcox_pvalue, "BH"), nrow = 3)
rownames(wilcox_adj_pvalue) <- rownames(wilcox_pvalue)
colnames(wilcox_adj_pvalue) <- colnames(wilcox_pvalue)
round(wilcox_adj_pvalue, 3)
round(wilcox_pvalue, 3)
### correct multipling testing for each contrast
wilcox_adj_pvalue_2 <- t(apply(wilcox_pvalue, 1, function(x) p.adjust(x, "BH")))
round(wilcox_adj_pvalue_2, 3)
# 3C - Pvalue plots
wilcox_pvalue_v <- reshape2::melt(wilcox_pvalue) %>%
    set_names(c("TumorType", "Type", "P_value")) %>%
    mutate(
        Sub = str_sub(Type, 3, 5),
        Type2 = str_c(str_sub(Type, 1, 1), "p", str_sub(Type, 7, 7))
    )
p$c <- ggplot(wilcox_pvalue_v, aes(Type2, -log10(P_value), fill = TumorType)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = .7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) +
    geom_hline(yintercept = -log10(0.05)) +
    labs(x = "Type", fill = "") +
    facet_grid(~Sub, space = "free", scales = "free") +
    theme(legend.position = "bottom")
p$c
ggsave(paste0(foldername, "/[3C TOP] WilcoxPvalue.pdf"), width = 11.5, height = 5)
wilcox_pvalue %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    write_tsv(paste0(foldername, "/[3C TOP] WilcoxPvalue.tsv"))
wilcox_adj_pvalue_v2 <- reshape2::melt(wilcox_adj_pvalue_2) %>%
    set_names(c("TumorType", "Type", "P_value")) %>%
    mutate(
        Sub = str_sub(Type, 3, 5),
        Type2 = str_c(str_sub(Type, 1, 1), "p", str_sub(Type, 7, 7))
    )
p_c2 <- ggplot(wilcox_adj_pvalue_v2, aes(Type2, -log10(P_value), fill = TumorType)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = .7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) +
    geom_hline(yintercept = -log10(0.1)) +
    labs(x = "Type", y = "-log10(FDR)", fill = "") +
    facet_grid(~Sub, space = "free", scales = "free") +
    theme(legend.position = "bottom")
p_c2
ggsave(paste0(foldername, "/[3C TOP] WilcoxAdjustedPvaluePerContrast.pdf"), width = 11.5, height = 5)
wilcox_adj_pvalue_2 %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    write_tsv(paste0(foldername, "/[3C TOP] WilcoxAdjustedPvaluePerContrast.tsv"))
# Bootstraping signatures
sanger_file <- "./signatures_probabilities.txt"
bootstrapping_rdata <- pathify("map19-bootstrappingResults.RData")
source("./calculateMutationSignature.R")
source("./getNearbyBases.R")
snv_nb_sam_h <- snv_nb_sam_h %>%
    rename(Tumor_Sample_Barcode = sam_id) %>%
    as.data.frame()
swt <- calGSW(sanger_file, snv_nb_sam_h)
dim(swt)
table(apply(swt, 1, sum))
summary(as.vector(swt))
swt <- swt[this_map_sam$sample_id, ]
if (!file.exists(bootstrapping_rdata)) {
    suppressWarnings({
        ## Perform boostrapping analysis
        sanger <- read.delim(sanger_file, header = TRUE, stringsAsFactors = FALSE)
        sanger <- sanger[1:96, 1:33]
        rownames(sanger) <- as.character(sanger$Somatic.Mutation.Type)
        # sanger[1:3, 1:5]
        sanger_mat <- as.matrix(sanger[, 4:33])
        rownames(sanger_mat) <- as.character(sanger$Somatic.Mutation.Type)
        sanger_mat[1:3, 1:3]

        snv_nb_mat <- t(as.matrix(snv_nb_sam_h[, 2:97]))
        colnames(snv_nb_mat) <- snv_nb_sam_h$Tumor_Sample_Barcode
        rownames(snv_nb_mat) <- gsub("\\(", "\\[", rownames(snv_nb_mat))
        rownames(snv_nb_mat) <- gsub("\\)", "\\]", rownames(snv_nb_mat))
        snv_nb_mat[1:4, 1:4]

        ## We matched the order of SNV between two matrixes.

        sanger_mat <- sanger_mat[rownames(snv_nb_mat), ]
        sanger <- sanger[rownames(snv_nb_mat), ]

        set.seed(123)
        boot_sanger_wt <- bootstrapGeneSignatureWeight(snv_nb_mat, sanger_mat, B = 2000)
        save(boot_sanger_wt, file = bootstrapping_rdata)
    })
} else {
    load(bootstrapping_rdata)
}
boot_sanger_wt_quantileTop5Percent <- sapply(
    boot_sanger_wt,
    function(x) apply(x, 2, quantile, prob = 0.95, na.rm = T)
)
class(boot_sanger_wt_quantileTop5Percent)
dim(boot_sanger_wt_quantileTop5Percent)
rownames(boot_sanger_wt_quantileTop5Percent) <- sub("ignature.", "", rownames(boot_sanger_wt_quantileTop5Percent))
boot_sanger_wt_quantileTop5Percent <- boot_sanger_wt_quantileTop5Percent[, rownames(swt)]

identical(rownames(swt), rownames(t(boot_sanger_wt_quantileTop5Percent)))
identical(colnames(swt), colnames(t(boot_sanger_wt_quantileTop5Percent)))

sanger_wt_sig <- swt > t(boot_sanger_wt_quantileTop5Percent)
dim(sanger_wt_sig)

# 3e - Heatmap of signature weights
ht_col <- circlize::colorRamp2(seq(0, .3, length = 3), c("blue", "#EEEEEE", "red"))

mean_wt_sig <- apply(swt, 2, mean)
mean_wt_sig

ha <- HeatmapAnnotation(
    mean_bar = anno_barplot(mean_wt_sig,
        border = F,
        axis = T,
        gp = gpar(fill = "gray48", col = "gray48")
    ),
    annotation_height = unit(2, "cm")
)

ht_e <- Heatmap(swt,
    name = "ht1",
    col = ht_col,
    cluster_rows = F,
    cluster_columns = F,
    cluster_row_slices = FALSE,
    clustering_distance_rows = "pearson",
    clustering_method_rows = "ward.D2",
    clustering_distance_columns = "pearson",
    clustering_method_columns = "ward.D2",
    column_title = "",
    column_names_gp = gpar(fontsize = 5),
    row_names_gp = gpar(fontsize = 5),
    row_title_gp = gpar(fontsize = 8),
    row_split = this_map_sam$type,
    heatmap_legend_param = list(color_bar = "continuous", title = "")
)

pdf(paste0(foldername, "/MSweightsHeatmapClustered.pdf"), width = 8.5, height = 3.5)
draw(ht_e, show_heatmap_legend = T)
dev.off()
saveRDS(ht_e, pathify("ht_e_clustered.Rds"))

draw(ht_e)

# 3d - barplot on the top of heatmap**
mean_wt_sig <- t(apply(swt, 2, function(x) tapply(x, this_map_sam$type, mean)))
mean_wt_sig

p$d <- as.data.frame(mean_wt_sig) %>%
    mutate(sid = factor(rownames(mean_wt_sig), levels = paste0("S", 1:30))) %>%
    gather(-sid, key = "Type", value = "Ave.Weight") %>%
    mutate(Type = factor(Type, levels = levels(this_map_sam$type))) %>%
    ggplot(mapping = aes(x = sid, y = Ave.Weight)) +
    geom_bar(aes(fill = Type),
        stat = "identity",
        position = position_dodge(width = 0.9),
        width = .7
    ) +
    scale_fill_manual(
        values = c("#3dbce7", "#ffac49", "black"),
        breaks = c("NU/LGIN", "HGIN", "UC")
    ) +
    xlab("") +
    ylab("Average Weight") +
    theme(legend.position = "top")
p$d

ggsave(paste0(foldername, "/MSweightsBarplot.pdf"), width = 8.5, height = 3)

# 3g - Heatmap of boostrapping result**
ht_col <- c("gray95", "steelblue")
names(ht_col) <- c(0, 1)

identical(rownames(sanger_wt_sig), rownames(swt))

ht_g <- Heatmap(sanger_wt_sig + 0,
    name = "ht-1",
    show_row_names = F, show_column_names = T,
    cluster_rows = F, cluster_columns = F,
    col = ht_col,
    row_order = unlist(row_order(ht_e)),
    clustering_method_rows = "ward.D", clustering_distance_rows = "euclidean", # "pearson",
    column_title = "",
    row_names_gp = gpar(fontsize = 3),
    row_title_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 5),
    split = this_map_sam$type,
    gap = unit(1, "mm"),
    heatmap_legend_param = list(color_bar = "discrete", title = "", legend_direction = "horizontal", nrow = 1)
)

row_order(ht_g)
row_order(ht_e)

pdf(paste0(foldername, "/SignificantGeneSignatureScore.pdf"), width = 8.5, height = 3.5)
draw(ht_g, show_heatmap_legend = T, heatmap_legend_side = "bottom")
dev.off()

draw(ht_g)

# Kruskal-Wallis test
kruskal.re <- apply(swt, 2, function(x) {
    tdat <- cbind(this_map_sam, weight = x)
    kruskal.test(weight ~ type, data = tdat)
})
kruskal.pvalue <- sapply(kruskal.re, function(x) x$p.value)

round(sort(kruskal.pvalue), 3)

adj.kruskal.pvalue <- p.adjust(kruskal.pvalue, method = "BH")
round(sort(adj.kruskal.pvalue), 3)

## 3f - Plot p-value
pdat <- data.frame(
    Signature = str_replace(names(kruskal.pvalue), "S", "") %>%
        parse_factor(),
    Pvalue = -log10(kruskal.pvalue)
)

p$f <- ggplot(pdat, aes(Signature, Pvalue)) +
    geom_bar(stat = "identity", width = .7) +
    geom_hline(aes(yintercept = -log10(.05)), colour = "#BB0000", linetype = "dashed") +
    ylab("-log10(P-value)")
p$f

ggsave(paste0(foldername, "/[3F] PvaluesforSignatures.pdf"), width = 11.5, height = 6)

pdat <- data.frame(
    Signature = str_replace(names(adj.kruskal.pvalue), "S", "") %>%
        parse_factor(),
    Pvalue = -log10(adj.kruskal.pvalue)
)

p_f2 <- ggplot(pdat, aes(Signature, Pvalue)) +
    geom_bar(stat = "identity", width = .7) +
    ylab("-log10(Adjusted P-value)") +
    geom_hline(aes(yintercept = -log10(.1)), colour = "#BB0000", linetype = "dashed")
p_f2

ggsave(paste0(foldername, "/[3F] AdjustedPvaluesforSignatures.pdf"), width = 11.5, height = 6)


## All plots
p$e <- grid.grabExpr(draw(ht_e))
p$g <- grid.grabExpr(draw(ht_g))
p <- p[c("a", "b", "c", "d", "e", "f", "g")]
layout <-
    "
  ABBBBDD
  ABBBBEE
  ABBBBFF
  ACCCCGG
"
wrap_plots(p, design = layout)
ggsave(pathify("map26.png"))
# Signatures by mutSignatures
## [MutSignatures: an R package for extraction and analysis of cancer mutational signatures](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7589488/)
## De novo extraction of signatures
# Count
blca_counts <- VAF_tbl %>%
    filter(VAF > treshold) %>%
    select(sam_id, mut_type) %>%
    as.data.frame() %>%
    countMutTypes(
        mutType_colName = "mut_type",
        sample_colName = "sam_id"
    )
blca_counts
# how many signatures should we extract?
num_sign <- 4
# Define parameters for the non-negative matrix factorization procedure.
# you should parallelize if possible
blca_params <-
    mutSignatures::setMutClusterParams(
        num_processesToExtract = num_sign, # num signatures to extract
        num_totIterations = 20, # bootstrapping: usually 500-1000
        num_parallelCores = 4
    ) # total num of cores to use (parallelization)

png(str_c("results/map19_mutSignatures-SilhouettePlot-n=", num_sign, "_t=", treshold, ".png"))
blca.analysis <- decipherMutationalProcesses(blca_counts, params = blca_params)
# Retrieve signatures (results)
map26_sigs <- blca.analysis$Results$signatures
# Retrieve exposures (results)
blca.exp <- blca.analysis$Results$exposures
# Plot signature 1 (standard barplot, you can pass extra args such as ylim)
m1 <- wrap_elements(~ msigPlot(map26_sigs, signature = 1, ylim = c(0, 0.10)))
m2 <- wrap_elements(~ msigPlot(map26_sigs, signature = 2, ylim = c(0, 0.10)))
m3 <- wrap_elements(~ msigPlot(map26_sigs, signature = 3, ylim = c(0, 0.10)))
m4 <- wrap_elements(~ msigPlot(map26_sigs, signature = 4, ylim = c(0, 0.10)))
### Might throw an error
wrap_plots(list(m1, m2, m3, m4), ncol = ggsave(str_c("results/map19_mutSignatures-signatures-n=", num_sign, "_t=", treshold, ".png")))
# Plot exposures (ggplot2 object, you can customize as any other ggplot2 object)
msigPlot(blca.exp, main = "BLCA samples") +
    scale_fill_manual(values = c("#1f78b4", "#cab2d6", "#ff7f00", "#a6cee3"))
ggsave(paste0(foldername, "/mutSigSamplesBarplot-n=", num_sign, ".png"))
dt <- blca.exp %>%
    {
        dt <- .@exposures
        colnames(dt) <- getSampleIdentifiers(.)
        rownames(dt) <- getSignatureIdentifiers(.)
        dt[, ur_samples$sample_id]
    }
mat <- t(scale(dt))
sample_classification_data <- samples_data %>%
    select(sample_id, group) %>%
    deframe() %>%
    .[rownames(mat)]
right_anno <- rowAnnotation(
    group = sample_classification_data,
    col = list(
        group = c(`NU/LGIN` = "#3dbce7", HGIN = "#ffac49", UC = "black")
    )
)
ht <- Heatmap(
    mat,
    # cluster_rows = FALSE,
    cluster_columns = FALSE,
    # row_split = ur_samples$group,
    right_annotation = right_anno
)
png(paste0(foldername, "/mutSigSamplesHeatmap-n=", num_sign, ".png"))
draw(ht)
dev.off()

draw(ht)

COSMIC_sigs <- mutSignatures::getCosmicSignatures()

COSMIC_sigs_2 <- mutSigData$blcaSIGS %>%
    dplyr::select(starts_with("COSMIC")) %>%
    as.mutation.signatures()

BLCA_knwn_sigs <- mutSigData$blcaSIGS %>%
    dplyr::select(starts_with("BLCA")) %>%
    as.mutation.signatures()

# Compare de-novo signatures with selected COSMIC signatures
msig_cosmic <- matchSignatures(
    mutSign = map26_sigs, reference = COSMIC_sigs,
    threshold = 0.45, plot = TRUE
)
msig_blca <- matchSignatures(
    mutSign = map26_sigs, reference = BLCA_knwn_sigs,
    threshold = 0.45, plot = TRUE
)

m1 <- msig_cosmic$plot + ggtitle("Match to COSMIC signs.")
m2 <- msig_blca$plot + ggtitle("Match to known BLCA signs.")
m1
ggsave(paste0(foldername, "/mutSignatures-ourSignatures-to-COSMIC-n=", num_sign, ".png"))
m2
ggsave(paste0(foldername, "/mutSignatures-ourSignatures-to-knownBLCA-n=", num_sign, ".png"))

## Estimate the contribution of a panel of known mutational signatures

# Run analysis
map19_COSMIC_sigs <- resolveMutSignatures(
    mutCountData = blca_counts,
    signFreqData = COSMIC_sigs
)

map19_BLCA_sigs <- resolveMutSignatures(
    mutCountData = blca_counts,
    signFreqData = BLCA_knwn_sigs
)

# Retrieve exposures (results)
blca.exp.1x <- map26_COSMIC_sigs$Results$count.result
blca.exp.2x <- map26_BLCA_sigs$Results$count.result

# Plot exposures
bp1 <- msigPlot(blca.exp.1x, main = "Map19 | COSMIC sigs.")
bp2 <- msigPlot(blca.exp.2x, main = "Map19 | pre-blca sigs.")

# Visualize
bp1 + bp2

dt <- map26_COSMIC_sigs$Results$count.result %>%
    {
        dt <- .@exposures
        colnames(dt) <- getSampleIdentifiers(.)
        rownames(dt) <- getSignatureIdentifiers(.)
        dt
    }

ht <- Heatmap(
    t(scale(dt))[ur_samples$sample_id, ],
    # t(dt)[ur_samples$sam_id, ],
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = ur_samples$group
)

png(paste0(foldername, "/mutSignatures-samples-COSMIC-heatmap-n=", num_sign, ".png"))
draw(ht)
dev.off()

draw(ht)

dt <- blca.exp %>%
    {
        dt <- .@exposures
        colnames(dt) <- getSampleIdentifiers(.)
        rownames(dt) <- getSignatureIdentifiers(.)
        dt[, ur_samples$sample_id]
    }
dt %>%
    # scale() %>%
    as.data.frame() %>%
    rownames_to_column("signature") %>%
    pivot_longer(-signature, names_to = "sample_id", values_to = "val") %>%
    mutate(
        row = str_match(sample_id, "[:alpha:]"),
        column = str_match(sample_id, "(?<=[:alpha:])[:digit:]+")
    ) %>%
    ggplot(aes(row, column, fill = val)) +
    geom_tile() +
    # scale_x_discrete(limits = rev) +
    # scale_y_discrete(limits = rev) +
    scale_fill_gradient(low = "white", high = "red") +
    facet_wrap(~signature)

blca.cosmic.exp <- map26_COSMIC_sigs$Results$count.result
dt <- blca.cosmic.exp %>%
    {
        dt <- .@exposures
        colnames(dt) <- getSampleIdentifiers(.)
        rownames(dt) <- getSignatureIdentifiers(.)
        dt[, ur_samples$sample_id]
    }
dt %>%
    rownames_to_column("signature") %>%
    filter(signature %in% str_c("COSMIC", c(1, 6, 12), sep = ".")) %>%
    pivot_longer(-signature, names_to = "sample_id", values_to = "val") %>%
    mutate(
        row = str_match(sample_id, "[:alpha:]"),
        column = str_match(sample_id, "(?<=[:alpha:])[:digit:]+")
    ) %>%
    ggplot(aes(row, column, fill = val)) +
    geom_tile() +
    # scale_x_discrete(limits = rev) +
    # scale_y_discrete(limits = rev) +
    scale_fill_gradient(low = "white", high = "red") +
    facet_wrap(~signature)
