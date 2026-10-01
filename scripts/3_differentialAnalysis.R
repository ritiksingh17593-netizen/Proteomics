# Differential expression analysis: dNLS vs WT
# Input:  data/processed/report.pg_matrix.tsv
# Outputs: differential_expression_results.csv
#          volcano_plot.png

# Load required packages
library(limma)
library(ggplot2)
library(ggrepel)
library(dplyr)

# Set the input and output file paths.
input_file <- "data/processed/report.pg_matrix.tsv"
results_file <- "differential_expression_results.csv"
plot_file <- "volcano_plot.png"

# Read the tab-separated protein group matrix.
# check.names = FALSE preserves the column names from the input file.
protein_data <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Check that the file has the expected metadata and sample columns.
if (ncol(protein_data) < 12) {
  stop("The input file must contain at least 12 columns.")
}
if (!"Protein.Group" %in% names(protein_data)) {
  stop("The input file must contain a 'Protein.Group' column.")
}

# Rename columns 7–12 to identify the three replicates in each group.
sample_names <- c("dNLS_1", "dNLS_2", "dNLS_3", "WT_1", "WT_2", "WT_3")
names(protein_data)[7:12] <- sample_names

# Preserve the protein IDs and metadata for joining to the statistical results.
protein_metadata <- protein_data %>%
  transmute(
    Protein.Group = as.character(Protein.Group),
    Genes = Genes,
    Description = First.Protein.Description
  )

# Extract the six sample columns and convert them to a numeric matrix.
intensity_matrix <- as.matrix(
  as.data.frame(
    lapply(protein_data[, sample_names, drop = FALSE], function(x) {
      as.numeric(as.character(x))
    }),
    check.names = FALSE
  )
)
colnames(intensity_matrix) <- sample_names
rownames(intensity_matrix) <- as.character(protein_data$Protein.Group)

# Treat zero intensity as missing, then log2-transform the measured intensities.
intensity_matrix[intensity_matrix == 0] <- NA
intensity_matrix <- log2(intensity_matrix)

# Retain proteins with at least two measured values in either experimental group.
dnls_samples <- c("dNLS_1", "dNLS_2", "dNLS_3")
wt_samples <- c("WT_1", "WT_2", "WT_3")

keep_proteins <-
  rowSums(!is.na(intensity_matrix[, dnls_samples, drop = FALSE])) >= 2 |
  rowSums(!is.na(intensity_matrix[, wt_samples, drop = FALSE])) >= 2

filtered_matrix <- intensity_matrix[keep_proteins, , drop = FALSE]

# Set up the group design. The contrast below estimates dNLS minus WT,
# so a positive logFC indicates higher abundance in dNLS.
group <- factor(
  c("dNLS", "dNLS", "dNLS", "WT", "WT", "WT"),
  levels = c("WT", "dNLS")
)

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# Fit a linear model for each protein, define the dNLS-versus-WT contrast,
# and apply empirical Bayes moderation to the standard errors.
fit <- lmFit(filtered_matrix, design)
contrast_matrix <- makeContrasts(dNLS - WT, levels = design)
fit_contrasted <- contrasts.fit(fit, contrast_matrix)
fit_ebayes <- eBayes(fit_contrasted)

# Extract results for all retained proteins, including log2 fold change,
# raw P value, and Benjamini–Hochberg adjusted P value.
limma_results <- topTable(
  fit_ebayes,
  coef = 1,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# Move protein IDs into a column and retain the requested result fields.
limma_results$Protein.Group <- rownames(limma_results)
limma_results <- limma_results %>%
  select(Protein.Group, logFC, P.Value, adj.P.Val)

# Join results to the protein metadata for a readable output table.
final_results <- limma_results %>%
  left_join(protein_metadata, by = "Protein.Group") %>%
  select(Protein.Group, Genes, Description, logFC, P.Value, adj.P.Val)

# Label proteins by significance using the requested FDR and fold-change cutoffs.
final_results <- final_results %>%
  mutate(
    Significance = case_when(
      !is.na(adj.P.Val) & adj.P.Val < 0.05 & logFC > 1  ~ "Upregulated",
      !is.na(adj.P.Val) & adj.P.Val < 0.05 & logFC < -1 ~ "Downregulated",
      TRUE                                               ~ "Not significant"
    )
  )

# Export the merged differential expression results.
write.csv(final_results, results_file, row.names = FALSE, na = "")

# Select up to ten significant proteins with the smallest adjusted P values
# for labeling. Use the protein ID if a gene name is unavailable.
label_data <- final_results %>%
  filter(Significance != "Not significant", !is.na(adj.P.Val)) %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 10) %>%
  mutate(
    Label = ifelse(is.na(Genes) | Genes == "", Protein.Group, Genes)
  )

# Prepare the y-axis values. The lower bound prevents log10(0) if an
# adjusted P value is reported as exactly zero.
plot_data <- final_results %>%
  filter(!is.na(logFC), !is.na(adj.P.Val)) %>%
  mutate(
    neg_log10_FDR = -log10(pmax(adj.P.Val, .Machine$double.xmin)),
    Significance = factor(
      Significance,
      levels = c("Downregulated", "Not significant", "Upregulated")
    )
  )

# Create the volcano plot.
volcano_plot <- ggplot(
  plot_data,
  aes(x = logFC, y = neg_log10_FDR, color = Significance)
) +
  geom_point(alpha = 0.75, size = 1.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_text_repel(
    data = label_data,
    aes(x = logFC, y = -log10(pmax(adj.P.Val, .Machine$double.xmin)),
        label = Label),
    inherit.aes = FALSE,
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c(
      "Downregulated" = "#3B82F6",
      "Not significant" = "grey70",
      "Upregulated" = "#DC2626"
    ),
    drop = FALSE
  ) +
  labs(
    title = "Differential expression: dNLS vs WT",
    x = expression(log[2] * " fold change (dNLS - WT)"),
    y = expression(-log[10] * " adjusted P value"),
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

# Save a high-resolution PNG of the volcano plot.
ggsave(
  filename = plot_file,
  plot = volcano_plot,
  width = 8,
  height = 7,
  units = "in",
  dpi = 600
)
