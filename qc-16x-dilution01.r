## ----globals -------------------------------------------------------------
packages <- c("eval.save", "metaboprep", "tidyverse", "knitr") 
lapply(packages, require, character.only=T)

dir <- paths
eval.save.dir(dir$cache)

## ----load.olink -------------------------------------------------------------

## read in the raw olink results so I can get LOD information
olink <- as.data.frame(data.table::fread(
			file.path(dir$data, 
			"GB390725-RB_1_16_Extended_2026-07-06.csv"))) 
colnames(olink) <- colnames(olink) |>
					make.names()|>
					tolower()
olink <- as_tibble(olink) |>
			filter(!grepl("CONTROL", sampleid))
olink <- olink |>
			mutate(patient.id     = sub("_.*", "", sampleid),
				dilution = ifelse(grepl("_", sampleid), sub(".*_", "", sampleid), 0)) |>
			mutate(dilution = factor(dilution, 
				levels = sort(unique(as.numeric(dilution)))))

## create dilution series lookup table
raw <- data.frame(sampleid = unique(olink$sampleid)) |>
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
unique(pheno.long$patient.id)

## ----dilution.ids -------------------------------------------------------------
pheno.long |>
	count(dilution) |>
	kable()

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

## ---- dist-other -------------------------------------------------------------
olink |>
	filter(dilution != "elisa_serum" & dilution != "elisa_pleural") |>
	filter(!assay %in% proteins) |>
	ggplot(aes(x = dilution, y = npx, fill = dilution)) +
	geom_violin(alpha = 0.5, colour = NA) +
	geom_boxplot(width = 0.15, outlier.size = 0.5, fill = "white", colour = "grey30") +
	facet_wrap(~ assay, scales = "free_y", ncol = 3) +
	labs(x = "Dilution", y = "NPX", fill = "Dilution") +
	theme_bw() +
	theme(
		legend.position = "bottom",
		axis.text.x     = element_text(angle = 45, hjust = 1),
		strip.text.x    = element_text(size = 7)
	)


## ---- il6-dilution-correlations -------------------------------------------------

# Calculate pairwise correlations between IL6 dilution levels
il6_data <- olink |>
  filter(assay == "IL6") 

dilution_levels <- levels(il6_data$dilution)
dil_pairs <- combn(dilution_levels, 2, simplify = FALSE)

cor_df <- bind_rows(lapply(dil_pairs, function(pair) {
  d1 <- pair[1]; d2 <- pair[2]

  # Individuals present at both dilutions
  patients_d1 <- il6_data |> filter(dilution == d1) |> pull(patient.id)
  patients_d2 <- il6_data |> filter(dilution == d2) |> pull(patient.id)
  common_patients <- intersect(patients_d1, patients_d2)

  # Get NPX values for common patients at each dilution
  df1 <- il6_data |>
    filter(dilution == d1, patient.id %in% common_patients) |>
    arrange(patient.id) |>
    select(patient.id, npx1 = npx)
  df2 <- il6_data |>
    filter(dilution == d2, patient.id %in% common_patients) |>
    arrange(patient.id) |>
    select(patient.id, npx2 = npx)

  pair_data <- inner_join(df1, df2, by = "patient.id")

  data.frame(
    dilution1 = d1,
    dilution2 = d2,
    pair      = paste(d1, "vs", d2),
    r         = cor(pair_data$npx1, pair_data$npx2, use = "pairwise.complete.obs"),
    n         = nrow(pair_data),
    stringsAsFactors = FALSE
  )
}))

cor_df |>
  arrange(desc(r)) |>
  kable(digits = 3, col.names = c("Dilution 1", "Dilution 2", "Pair", "Pearson r", "N individuals"))


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
dilution_scatter("IL6", olink = il6_data)

