## ----globals -------------------------------------------------------------
packages <- c("eval.save", "tidyverse") 
lapply(packages, require, character.only=T)

# set dirs  
dir <- paths
eval.save.dir(dir$cache)

## ----load.elisa -------------------------------------------------------------
## load raw clincal phenotype info 
elisa <-  data.table::fread(file.path(dir$data, 
                "il6_clean_los.txt")) 
colnames(elisa) <- colnames(elisa) |>
                make.names()|>
                tolower()
elisa <- as_tibble(elisa)				
str(elisa)

## ----check.elisa -------------------------------------------------------------
elisa|>
    summarise(
		n_elisa_serum = sum(!is.na(serum_il6)),
		n_elisa_pleural = sum(!is.na(pleural_il6)),
        n_both = sum(!is.na(serum_il6) & !is.na(pleural_il6), na.rm = T)
		) |>
		kable()
#	| n_elisa_serum| n_elisa_pleural| n_both|
#	|-------------:|---------------:|------:|
#	|            84|              86|     84|

## ----load.olink -------------------------------------------------------------
## get il6 feature_id - protein info
annot <- eval.ret(paste("annot", "pleural", sep=".")) |>
			select(feature_id, assay) |>
			mutate(protein = make.names(tolower(assay)))|>
			as.data.frame()

idx.il6 <- unlist(annot[which(annot$protein %in% "il6"),] )
#    feature_id assay protein
# 11   OID00482   IL6     il6

## proteins
olink.plasma <- 
	eval.ret(paste("prot.mat", "plasma", sep=".")) |>
		t() |>
		as.data.frame()|>
		rownames_to_column(var = "sample") |>
		as_tibble() |>
		select(sample, all_of(idx.il6["feature_id"]))|>
		rename(plasma = feature_id)|>
		mutate(olink.plasma = 1)
		
olink <- eval.ret(paste("prot.mat", "pleural", sep=".")) |>
		t() |>
		as.data.frame()|>
		rownames_to_column(var = "sample") |>
		as_tibble() |>
		select(sample, all_of(idx.il6["feature_id"]))|>
		rename(pleural = feature_id)|>
		mutate(olink.pleural = 1) |>
		left_join(olink.plasma, by = "sample")

## ----check.olink -------------------------------------------------------------
olink |>
    summarise(
		n_olink_pleural = sum(olink.pleural),
		n_olink_plasma = sum(olink.plasma, na.rm=T),
        n_both = sum(olink.pleural == 1 & olink.plasma == 1, na.rm = T)
		) |>
		kable()
#	| n_olink_pleural| n_olink_plasma| n_both|
#	|---------------:|--------------:|------:|
#	|             256|             80|     80|

## ----overlap -------------------------------------------------------------
addmargins(table(olink$sample %in% elisa$patient_id)) |> kable()

## ----join -------------------------------------------------------------
il6 <- olink |>
        full_join(elisa, by = c("sample" = "patient_id"))

## ----count -------------------------------------------------------------
il6 |>
    summarise(
		n_olink_pleural = sum(olink.pleural, na.rm=T),
		n_olink_plasma = sum(olink.plasma, na.rm=T),
        n_olink_both = sum(olink.pleural == 1 & olink.plasma == 1, na.rm = T),
		n_elisa_serum = sum(!is.na(serum_il6)),
		n_elisa_pleural = sum(!is.na(pleural_il6)),
		n_elisa_both = sum(!is.na(serum_il6) & !is.na(pleural_il6), na.rm = T)
		) |>
		kable()

