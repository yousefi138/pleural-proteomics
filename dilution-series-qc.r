## ----globals -------------------------------------------------------------
packages <- c("eval.save", "metaboprep", "tidyverse", "knitr") 
lapply(packages, require, character.only=T)

dir <- paths
eval.save.dir(dir$cache)

## ----load.data -------------------------------------------------------------

## load proetins
prot <- eval.ret(paste("prot.mat", "pleural-dilution-series", sep="."))

## ----access.pheno -------------------------------------------------------------

## create dilution series lookup table
raw <- data.frame(sampleid = colnames(prot)) |>
	mutate(
	patient.id     = sub("_.*", "", sampleid),
	dilution = ifelse(grepl("_", sampleid), sub(".*_", "", sampleid), 0)
	)

# restrict pheno to the n=16 samples included in the dilution series
pheno <- eval.ret("pheno") |>
			filter(patient.id %in% unique(raw$patient.id))

pheno.long <- raw |>
			left_join(pheno, by = "patient.id") |>
			mutate(dilution = factor(dilution, 
				levels = sort(unique(as.numeric(dilution)))))

## ----dilutions -------------------------------------------------------------
pheno.long |>
	count(dilution) |>
	kable()

## ----dilution.ids -------------------------------------------------------------
pheno.long |>
	split(pheno.long$dilution) |>
	map(~ .x$patient.id)

## ---- dilution-pair-correlations -------------------------------------------------------------

# All pairwise combinations of dilution levels
dilution_levels <- levels(pheno.long$dilution)
dil_pairs <- combn(dilution_levels, 2, simplify = FALSE)

cor_df <- bind_rows(lapply(dil_pairs, function(pair) {
  d1 <- pair[1]; d2 <- pair[2]

  # Individuals present at both dilutions
  patients_d1 <- pheno.long |> filter(dilution == d1) |> pull(patient.id)
  patients_d2 <- pheno.long |> filter(dilution == d2) |> pull(patient.id)
  common_patients <- intersect(patients_d1, patients_d2)

  # Sample IDs for common patients, ordered identically so columns align
  s1 <- pheno.long |>
    filter(dilution == d1, patient.id %in% common_patients) |>
    arrange(patient.id) |> pull(sampleid)
  s2 <- pheno.long |>
    filter(dilution == d2, patient.id %in% common_patients) |>
    arrange(patient.id) |> pull(sampleid)

  mat1 <- prot[, s1, drop = FALSE]
  mat2 <- prot[, s2, drop = FALSE]

  # Per-protein Pearson r across shared individuals
  r_vec <- vapply(rownames(prot), function(pid) {
    x <- as.numeric(mat1[pid, ])
    y <- as.numeric(mat2[pid, ])
    if (sum(complete.cases(x, y)) < 3L) return(NA_real_)
    cor(x, y, use = "pairwise.complete.obs", method = "pearson")
  }, numeric(1L))

  data.frame(
    dilution1 = d1,
    dilution2 = d2,
    pair      = paste(d1, "vs", d2),
    protein   = names(r_vec),
    r         = r_vec,
    n         = length(common_patients),
    stringsAsFactors = FALSE
  )
}))

## summary: median r per dilution pair
cor_summary <- cor_df |>
  group_by(pair, dilution1, dilution2, n) |>
  summarise(
    median_r = median(r, na.rm = TRUE),
    mean_r   = mean(r,   na.rm = TRUE),
    n_prot   = sum(!is.na(r)),
    .groups  = "drop"
  )
print(cor_summary)

## violin / box plot of per-protein r by dilution pair
p_dil_cor <- ggplot(cor_df, aes(x = pair, y = r)) +
  geom_violin(fill = "steelblue", alpha = 0.4, colour = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
  labs(
    title = "Per-protein correlation between dilution levels",
    x     = "Dilution pair",
    y     = "Pearson r (across individuals)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_dil_cor)

## ---- -------------------------------------------------------------

