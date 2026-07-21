## ----globals -------------------------------------------------------------
packages <- c("eval.save", "tidyverse", "knitr", "GGally") 
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

## ----load.elisa -------------------------------------------------------------
## load raw clincal phenotype info 
elisa <-  data.table::fread(file.path(dir$data, 
                "il6_clean_los.txt")) 
colnames(elisa) <- colnames(elisa) |>
                make.names()|>
                tolower()
elisa <- as_tibble(elisa) |>		
            select(patient_id, pleural_il6, serum_il6) |>
			mutate(log_pleural_il6 = log(pleural_il6))		
str(elisa)

## ----il6 -------------------------------------------------------------
il6 <- olink |>
        filter(olinkid == "OID00482") |>
        select(patient.id, npx) |>
        rename(olink_il6 = npx) |>
        full_join(elisa, by = c("patient.id" = "patient_id"))


## ----lod -------------------------------------------------------------
qc <- 
	olink|>
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

## ----count -------------------------------------------------------------
il6 |>
    summarise(
		n_olink_and_elisa = sum(!is.na(olink_il6) & !is.na(pleural_il6), na.rm=T)
		) |>
		kable()

## ----lm -------------------------------------------------------------
il6 |>
    select(olink_il6, pleural_il6, serum_il6, log_pleural_il6)|>
	ggpairs(lower = list(continuous = wrap("smooth", method = "lm")))
