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
# ===================================== HEATMAPS OF NON-SILENT MUTATIONS
# ======================================================================
suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(furrr)
    library(pins)
    library(tidyverse)
})
#---Input field data and group into categories (NU/LGIN, HGIN/UC)
samples_data <- read_tsv("samples_data.tsv") |>
    mutate(
        group2 = case_when(
            group4 %in% c("NU", "LGIN") ~ "NU/LGIN",
            group4 %in% c("HGIN", "UC") ~ "HGIN/UC",
            TRUE ~ NA_character_
        ) |>
            parse_factor(levels = c("NU/LGIN", "HGIN/UC")),
        group3 = parse_factor(group3, levels = c("NU/LGIN", "HGIN", "UC")),
        group4 = parse_factor(group4, levels = c("NU", "LGIN", "HGIN", "UC")),
        code = parse_factor(code, levels = c("N", "D1", "D2", "D3", "D4", "D5"))
    )
#---Input mutation data and filter for samples in consideration
vcf_data <- readRDS("vcf_data.final_unfiltered.rds") |>
    filter(sample_id %in% samples_data$sample_id)
#---Check unique mutation counts of each type in each map
result <- vcf_data |>
    select(patient_id, mutation_id, alt, ref, effect, variant_type) |>
    unique() |>
    group_by(patient_id, variant_type) |>
    count() |>
    pivot_wider(names_from = patient_id, values_from = n, values_fill = 0, id_expand = TRUE)
cat("Mutation counts in each map:\n")
print(result)
#---Plot nonsilent & null mutation counts in each map
plot_numbers_of_nonsilent_variants <- function(dt) {
    dt |>
        filter(effect %in% c("nonsilent", "null")) |>
        select(patient_id, mutation_id, effect) |>
        unique() |>
        ggplot(aes(patient_id, fill = effect)) +
        geom_bar()
}
pdf(paste0(foldername, "/Nonsilent variant counts.pdf"))
print(plot_numbers_of_nonsilent_variants(vcf_data))
dev.off()
#---Keep only mutations with mean normal VAF < 0.01
vcf_data_filtered <- vcf_data |>
    filter(is.na(mean_normal_VAF) | mean_normal_VAF < 0.01)
pdf(paste0(foldername, "/Nonsilent variant counts after filtering.pdf"))
print(plot_numbers_of_nonsilent_variants(vcf_data_filtered))
dev.off()
#---Keep only mutations with mean normal DP > 50
vcf_data_filtered_2 <- vcf_data_filtered |>
    filter(is.na(mean_normal_DP) | mean_normal_DP > 50)
pdf(paste0(foldername, "/Nonsilent variant counts after filtering 2.pdf"))
print(plot_numbers_of_nonsilent_variants(vcf_data_filtered_2))
dev.off()
#---Plot sequencing depths of mutations in each map
plot_sequencing_depths <- function(dt) {
    dt |>
        ggplot(aes(DP, color = sample_id)) +
        geom_density(alpha = 1) +
        scale_x_log10() +
        facet_wrap(~patient_id) +
        theme(legend.position = "none")
}
#---Calculate thresholds on mean sequencing depths for each sample:
#---Threshold=max(20, mean log10(DP) - 2 sd for coding-region variants)
#---Variants are kept if they pass the threshold in at least 1 sample
per_sample_thresholds <- vcf_data_filtered_2 |>
    filter(effect != "noncoding") |>
    mutate(logDP = log10(DP)) |>
    group_by(patient_id, sample_id) |>
    summarize(
        DP_mean = mean(DP),
        DP_median = median(DP),
        logDP_mean = mean(logDP),
        logDP_sd = sd(logDP)
    ) |>
    mutate(
        threshold = 10^(logDP_mean - 2 * logDP_sd),
        threshold = pmax(threshold, 20)
    ) |>
    left_join(samples_data |> select(sample_id, batch))
pdf(paste0(foldername, "/Sequencing depth after filtering 2.pdf"))
print(vcf_data_filtered_2 |>
    filter(effect != "noncoding") |>
    plot_sequencing_depths() +
    geom_vline(aes(xintercept = threshold), data = per_sample_thresholds, linewidth = 0.3, alpha = 0.5) +
    geom_vline(xintercept = 50, linewidth = 0.3, linetype = "dashed"))
dev.off()
keep <- vcf_data_filtered_2 |>
    left_join(per_sample_thresholds) |>
    filter(DP > threshold) |>
    select(patient_id, mutation_id) |>
    unique()
vcf_data_filtered_3 <- vcf_data_filtered_2 |>
    semi_join(keep)
filtered_nonsilent_snvs <- vcf_data_filtered_3 |>
    filter(effect %in% c("nonsilent", "null"))
pdf(paste0(foldername, "/Sequencing depth after filtering 3.pdf"))
print(plot_sequencing_depths(filtered_nonsilent_snvs))
dev.off()
pdf(paste0(foldername, "/Sequencing depth after filtering 3.pdf"))
print(plot_sequencing_depths(filtered_nonsilent_snvs))
dev.off()
pdf(paste0(foldername, "/Nonsilent variant counts after filtering 3.pdf"))
print(plot_numbers_of_nonsilent_variants(filtered_nonsilent_snvs))
dev.off()
vcf_data_filtered_3 |>
    select(patient_id, mutation_id, alt, ref, effect, variant_type) |>
    unique() |>
    group_by(patient_id, variant_type) |>
    count() |>
    pivot_wider(names_from = patient_id, values_from = n, values_fill = 0, id_expand = TRUE)
saveRDS(vcf_data_filtered_3, paste0(foldername, "/filtered_variants.Rds"))
saveRDS(filtered_nonsilent_snvs, paste0(foldername, "/filtered_nonsilent_variants.Rds"))
filtered_variants <- readRDS(paste0(foldername, "/filtered_variants.Rds"))
filtered_nonsilent_variants <- readRDS(paste0(foldername, "/filtered_nonsilent_variants.Rds"))
#---Classify each mutation by number of samples it spreads across
dt <- filtered_nonsilent_variants |>
    group_by(patient_id, mutation_id) |>
    mutate(
        n = n(),
        spread = case_when(
            n == 1 ~ "private",
            n > 1 & n <= 10 ~ "2 - 10 samples",
            n > 10 & n <= 20 ~ "11 - 20 samples",
            n > 20 & n <= 30 ~ "21 - 30 samples",
            n > 30 ~ "31+ samples",
            TRUE ~ "ERROR"
        )
    ) |>
    left_join(samples_data)
dt <- dt |>
    ungroup() |>
    mutate(
        spread = parse_factor(spread, levels = names(colors$variant_spread))
    )
#---Plot the SFS in each field for each map
x <- dt |>
    nest_by(patient_id) |>
    deframe()
plot <- function(dt, name) {
    bkg <- dt |>
        select(sample_id, row, column, group4) |>
        unique()
    dt |>
        ggplot() +
        geom_histogram(aes(VAF, fill = spread)) +
        geom_point(aes(color = group4), x = Inf, y = Inf, data = bkg, size = 10) +
        facet_grid(row ~ column) +
        scale_fill_manual(values = colors$variant_spread) +
        scale_color_manual(values = colors$blca_4_class_red) +
        labs(title = name)
}
pdf(paste0(foldername, "/SFS_by_field.pdf"), width = 15, height = 10)
print(x |> imap(plot))
dev.off()
#---Prepare materials for plotting heatmaps
heatmaps_data <- filtered_nonsilent_variants |>
    select(patient_id, sample_id, mutation_id, VAF) |>
    nest_by(patient_id) |>
    deframe()
mats <- heatmaps_data |>
    map(pivot_wider, names_from = sample_id, values_from = VAF, values_fill = 0) |>
    map(column_to_rownames, var = "mutation_id") |>
    map(as.matrix)
plot_heatmap <- function(mat,
                         samples_data,
                         name,
                         ...,
                         split_columns = FALSE,
                         col_annotations = list(
                             Group = colors$blca_4_class_red,
                             Grade = colors$blca_codes
                         ),
                         max_rows = 10000) {
    if (nrow(mat) > max_rows) {
        message("Number of rows in matrix ", name, " exceeds max_rows:", max_rows, ". Subsampling...")
        keep <- sample(1:nrow(mat), size = max_rows)
        mat <- mat[keep, ]
    }
    colors1 <- circlize::colorRamp2(
        seq(0, 0.5, length = 3),
        c("#3952a4", "#eff0f1", "#ee1f25")
    )
    column_split <- if (is.character(split_columns)) {
        x <- samples_data |>
            select(sample_id, !!sym(split_columns)) |>
            deframe()
        x[colnames(mat)]
    } else {
        NULL
    }
    Group <- samples_data %>%
        select(sample_id, group4) %>%
        deframe() %>%
        .[colnames(mat)]
    Classification <- samples_data %>%
        select(sample_id, code) %>%
        deframe() %>%
        .[colnames(mat)]
    samples_anno <- columnAnnotation(
        Group = Group,
        Grade = Classification,
        col = list(
            Group = col_annotations$Group,
            Grade = col_annotations$Grade
        )
    )
    row_anno <- rowAnnotation(n = anno_barplot(rowSums(mat > 0.01)))
    ht <- Heatmap(
        mat,
        name = name,
        show_row_names = FALSE,
        column_names_gp = gpar(fontsize = 10),
        col = colors1,
        use_raster = TRUE,
        top_annotation = samples_anno,
        right_annotation = row_anno,
        column_split = column_split,
        ...
    )
    draw(ht)
}
#---Plot VAF heatmaps
ht1 <- mats %>%
    imap(~ plot_heatmap(
        .x,
        samples_data,
        name = .y,
        clustering_method_rows = "ward.D"
    ))
pdf(paste0(foldername, "/Heatmaps - raw.pdf"), width = 8, height = 8)
walk(ht1, print)
dev.off()
ht1b <- mats %>%
    imap(~ plot_heatmap(
        .x,
        samples_data,
        name = .y,
        clustering_method_rows = "ward.D",
        split_columns = "group2",
        cluster_column_slices = FALSE
    ))
pdf(paste0(foldername, "/[2A BOTTOM] Heatmaps - raw & splitted.pdf"), width = 8, height = 8)
walk(ht1b, draw)
dev.off()
#---Keep only variants present in more than 1 sample with VAF > 0.01
keep <- filtered_nonsilent_variants |>
    group_by(patient_id, mutation_id) |>
    filter(sum(VAF > 0.01) > 1) |>
    select(patient_id, mutation_id) |>
    unique()
heatmaps_data2 <- filtered_nonsilent_variants |>
    semi_join(keep) |>
    select(patient_id, sample_id, mutation_id, VAF) |>
    nest_by(patient_id) |>
    deframe()
mats2 <- heatmaps_data2 |>
    map(pivot_wider, names_from = sample_id, values_from = VAF, values_fill = 0) |>
    map(column_to_rownames, var = "mutation_id") |>
    map(as.matrix)
ht2 <- mats2 %>%
    imap(~ plot_heatmap(
        .x,
        samples_data,
        name = .y,
        row_split = 2,
        clustering_method_rows = "ward.D",
        clustering_distance_rows = if (.y == "M21") "canberra" else "euclidean"
    ))
pdf(paste0(foldername, "/Heatmaps - filtered.pdf"), width = 8, height = 8)
walk(ht2, draw)
dev.off()
ht2b <- mats2 %>%
    imap(~ plot_heatmap(
        .x,
        samples_data,
        name = .y,
        row_split = 2,
        clustering_method_rows = "ward.D",
        clustering_distance_rows = if (.y == "M21") "canberra" else "euclidean",
        split_columns = "group2",
        cluster_column_slices = FALSE
    ))
pdf(paste0(foldername, "/[2B BOTTOM] Heatmaps - filtered & splitted.pdf"), width = 8, height = 8)
walk(ht2b, draw)
dev.off()
#---Plot mutation count by category for each map
pdf(paste0(foldername, "/Mutation count by category.pdf"))
print(filtered_nonsilent_variants |>
    group_by(patient_id, sample_id) |>
    count() |>
    left_join(samples_data) |>
    ggplot(aes(group3, n, color = group3, fill = group3)) +
    geom_boxplot(alpha = 0.5) +
    scale_color_manual(values = colors$blca_3_class_orig) +
    scale_fill_manual(values = colors$blca_3_class_orig) +
    facet_wrap(~patient_id, scales = "free", nrow = 2) +
    labs(
        y = "Number of mutations per sample",
        color = "Sample\nclassification",
        fill = "Sample\nclassification"
    ) +
    theme(axis.text.x = element_blank()))
dev.off()
# #---Save results
saveRDS(ht1, paste0(foldername, "/ht1.Rds"))
saveRDS(ht1b, paste0(foldername, "/ht1_splitted.Rds"))
saveRDS(ht2, paste0(foldername, "/ht2.Rds"))
saveRDS(ht2b, paste0(foldername, "/ht2_splitted.Rds"))
