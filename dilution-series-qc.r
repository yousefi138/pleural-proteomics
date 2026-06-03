## ----globals -------------------------------------------------------------
packages <- c("eval.save", "metaboprep", "tidyverse", "knitr") 
lapply(packages, require, character.only=T)

dir <- paths
eval.save.dir(dir$cache)

## ----load.olink -------------------------------------------------------------

## load proetins
prot <- eval.ret(paste("prot.mat", "pleural-dilution-series", sep="."))

## read in the raw olink results so I can get LOD information
olink <- as.data.frame(data.table::fread(
			file.path(dir$data, 
			"GB390725-RB_dilutionseries_Extended_2026-05-07.csv"))) 
colnames(olink) <- colnames(olink) |>
					make.names()|>
					tolower()
olink <- as_tibble(olink)
olink <- olink |>
			mutate(patient.id     = sub("_.*", "", sampleid),
				dilution = ifelse(grepl("_", sampleid), sub(".*_", "", sampleid), 0)) |>
			mutate(dilution = factor(dilution, 
				levels = sort(unique(as.numeric(dilution)))))

## create dilution series lookup table
raw <- data.frame(sampleid = colnames(prot)) |>
	mutate(
	patient.id     = sub("_.*", "", sampleid),
	dilution = ifelse(grepl("_", sampleid), sub(".*_", "", sampleid), 0)
	)

## ----load.elisa -------------------------------------------------------------
## load raw clincal phenotype info 
elisa <-  data.table::fread(file.path(dir$data, 
                "il6_clean_los.txt")) 
colnames(elisa) <- colnames(elisa) |>
                make.names()|>
                tolower()
elisa <- as_tibble(elisa) |>
			select(patient_id, serum_il6, pleural_il6)|>
			rename(patient.id = patient_id,
				elisa_serum = serum_il6,
				elisa_pleural = pleural_il6) |>
			filter(patient.id %in% unique(raw$patient.id)) |>
			pivot_longer(
				cols         = c(elisa_serum, elisa_pleural),
				names_to     = "dilution",
				values_to    = "npx"
			) |>
			mutate(dilution = factor(dilution),
					assay = "IL6")

# add to olink data
olink <- bind_rows(olink, elisa)

## ----access.pheno -------------------------------------------------------------

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

## ----lod -------------------------------------------------------------
qc <- 
	olink|>
		filter(dilution!="elisa_serum"&dilution!="elisa_pleural") |>
		group_by(assay, dilution)|>
		summarise(
			n = n(),
			n.belowlod = sum(belowlod, na.rm = T),
			freq.belowlod = {
				sum(belowlod, na.rm = T)/n() 
			},
			n.miss = sum(is.na(npx)),
			freq.miss = sum(is.na(npx))/n()
		)	

qc |>
	ggplot(aes(y = reorder(assay, n.belowlod), x = n.belowlod)) +
	geom_col() +
	facet_wrap(~ dilution) +
	labs(y = "Protein") 

## ----lod.sub -------------------------------------------------------------
proteins <- c("IL6","TNF", "IL8", "uPA", "CXCL1", "IFN-gamma", "CXCL5", "IL-17A", "IL-1 alpha", "CCL20")

qc <- 
	olink|>
		filter(dilution!="elisa_serum"&dilution!="elisa_pleural") |>
		filter(assay %in% proteins) |>
		group_by(assay, dilution)|>
		summarise(
			n = n(),
			n.belowlod = sum(belowlod, na.rm = T),
			freq.belowlod = {
				sum(belowlod, na.rm = T)/n() 
			},
			n.miss = sum(is.na(npx)),
			freq.miss = sum(is.na(npx))/n()
		)	

qc |>
	ggplot(aes(y = reorder(assay, n.belowlod), x = n.belowlod)) +
	geom_col() +
	facet_wrap(~ dilution) +
	labs(y = "Protein") 

## ---- dist-subset -------------------------------------------------------------
olink |>
	filter(dilution != "elisa_serum" & dilution != "elisa_pleural") |>
	filter(assay %in% proteins) |>
	ggplot(aes(x = dilution, y = npx, fill = dilution)) +
	geom_violin(alpha = 0.5, colour = NA) +
	geom_boxplot(width = 0.15, outlier.size = 0.5, fill = "white", colour = "grey30") +
	facet_wrap(~ assay, scales = "free_y", ncol =3) +
	labs(x = "Dilution", y = "NPX", fill = "Dilution") +
	theme_bw() +
	theme(
		legend.position = "bottom",
		axis.text.x     = element_text(angle = 45, hjust = 1),
		strip.text.x    = element_text(size = 7)
	)

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

# Summary: median r per dilution pair
cor_summary <- cor_df |>
  group_by(pair, dilution1, dilution2, n) |>
  summarise(
    median_r = median(r, na.rm = TRUE),
    mean_r   = mean(r,   na.rm = TRUE),
    n_prot   = sum(!is.na(r)),
    .groups  = "drop"
  )

cor_summary |>
  arrange(desc(mean_r)) |>
  kable(digits = 3, col.names = c("Pair", "Dilution 1", "Dilution 2", "N individuals", "Median r", "Mean r", "N proteins"))

## ---- dilution-pair-correlations-plot -------------------------------------------------------------

pair_order <- cor_summary |> arrange(desc(median_r)) |> pull(pair)
cor_df <- cor_df |> mutate(pair = factor(pair, levels = pair_order))

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

p_dil_cor

## ---- dilution-scatter-fn -------------------------------------------------------------

dilution_scatter <- function(protein, olink) {
	prot_data <- olink |> filter(assay == protein) |> droplevels()
	dil_levels <- levels(prot_data$dilution)
	dil_pairs  <- combn(dil_levels, 2, simplify = FALSE)

	prot_pairs <- bind_rows(lapply(dil_pairs, function(pair) {
		d1 <- pair[1]; d2 <- pair[2]
		df1 <- prot_data |> filter(dilution == d1) |> select(patient.id, npx1 = npx)
		df2 <- prot_data |> filter(dilution == d2) |> select(patient.id, npx2 = npx)
		inner_join(df1, df2, by = "patient.id") |>
			mutate(
				pair      = paste("Dilution", d1, "vs", d2),
				dilution1 = d1,
				dilution2 = d2
			)
	}))

	cor_labels <- prot_pairs |>
		group_by(pair, dilution1, dilution2) |>
		summarise(
			r     = cor(npx1, npx2, use = "pairwise.complete.obs"),
			label = paste0("r = ", round(r, 3)),
			.groups = "drop"
		)

	plots <- lapply(seq_along(dil_pairs), function(i) {
		p       <- dil_pairs[[i]]
		d1      <- p[1]; d2 <- p[2]
		lbl     <- paste("Dilution", d1, "vs", d2)
		df      <- prot_pairs  |> filter(pair == lbl)
		r_label <- cor_labels  |> filter(pair == lbl) |> pull(label)

		ggplot(df, aes(x = npx1, y = npx2)) +
			geom_point(alpha = 0.7) +
			geom_smooth(aes(colour = "Linear fit"), method = "lm", se = FALSE) +
			geom_abline(aes(intercept = 0, slope = 1, colour = "Identity (y = x)"),
				linetype = "dashed") +
			scale_colour_manual(
				name   = NULL,
				values = c("Linear fit" = "steelblue", "Identity (y = x)" = "grey40")
			) +
			guides(colour = guide_legend(
				override.aes = list(linetype = c("dashed", "solid"))
			)) +
			annotate("text", x = -Inf, y = Inf, label = r_label,
				hjust = -0.1, vjust = 1.4, size = 3.5) +
			labs(
				title = lbl,
				x     = paste0("NPX (Dilution ", d1, ")"),
				y     = paste0("NPX (Dilution ", d2, ")")
			) +
			theme_bw()
	})

	patchwork::wrap_plots(plots, ncol = 3, guides = "collect") +
		patchwork::plot_annotation(title = protein) &
		theme(legend.position = "bottom")
}

## ---- dilution-scatter-il6 -------------------------------------------------------------
dilution_scatter("IL6", olink = olink)

## ---- dilution-scatter-plots -------------------------------------------------------------
proteins <- c("TNF", "IL8", "uPA", "CXCL1", "IFN-gamma", "CXCL5", "IL-17A", "IL-1 alpha", "CCL20")
lapply(proteins, dilution_scatter, olink = olink)

## ---- linearity-check -------------------------------------------------------------
proteins10 <- c("IL6", "TNF", "IL8", "uPA", "CXCL1", "IFN-gamma", "CXCL5", "IL-17A", "IL-1 alpha", "CCL20")

# Adjacent dilution pairs and their expected NPX shift (log2 of fold-ratio)
# Undiluted is encoded as "0" in the data but represents a 1x factor
adj_pairs <- list(
	c("0",   "4"),
	c("4",   "16"),
	c("16",  "64"),
	c("64",  "256"),
	c("256", "1025")
)

dil_factor <- function(d) ifelse(d == "0", 1, as.numeric(d))

linearity_df <- bind_rows(lapply(adj_pairs, function(pair) {
	d1 <- pair[1]; d2 <- pair[2]
	expected <- log2(dil_factor(d2) / dil_factor(d1))

	df1 <- olink |>
		filter(dilution == d1, assay %in% proteins10) |>
		select(patient.id, assay, npx_d1 = npx)
	df2 <- olink |>
		filter(dilution == d2, assay %in% proteins10) |>
		select(patient.id, assay, npx_d2 = npx)

	inner_join(df1, df2, by = c("patient.id", "assay")) |>
		mutate(
			pair           = paste0(d1, "\u00d7 \u2192 ", d2, "\u00d7"),
			observed_shift = npx_d1 - npx_d2,
			expected_shift = expected
		)
}))

# Ordered factor so panels appear in dilution sequence
linearity_df <- linearity_df |>
	mutate(pair = factor(pair, levels = unique(pair)))

ggplot(linearity_df, aes(x = pair, y = observed_shift)) +
	geom_hline(aes(yintercept = expected_shift), 
			   linetype = "dashed", colour = "firebrick", linewidth = 0.7) +
	geom_boxplot(fill = "steelblue", alpha = 0.4, outlier.size = 0.8) +
	geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
	facet_wrap(~ assay, ncol = 2) +
	labs(
		x = "Adjacent dilution pair",
		y = "Observed NPX shift (lower \u2212 higher dilution)",
		caption = "Red dashed line = expected shift under perfect linearity (log2 of fold-ratio = 2 for all 4\u00d7 steps)"
	) +
	theme_bw() +
	theme(
		axis.text.x  = element_text(angle = 45, hjust = 1),
		strip.text   = element_text(size = 8),
		plot.caption = element_text(size = 7, colour = "grey40")
	)



