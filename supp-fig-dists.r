## ----fig-assoc -------------------------------------------------------------
packages <- c("eval.save", "dplyr", "tidyr", "ggplot2") 
lapply(packages, require, character.only=T)

dir <- paths
eval.save.dir(dir$cache)

## ----access.pheno -------------------------------------------------------------

# Add batch data and restrict pheno to those with proteins passing qc
pheno <- eval.ret("pheno") |>
			left_join(
				eval.ret(paste("batch", project, sep=".")), 
					by = c("patient.id" = "sample_id")) |>  

           ## keep only samples passing qc
            filter(excluded == FALSE) |>
			filter(sampleqc == "PASS")

## ----get protins -------------------------------------------------------------
prot.mat <- eval.ret(paste("prot.mat", project, sep="."))
prot <- prot.mat[,match(pheno$patient.id, colnames(prot.mat))]

annot <- eval.ret(paste("annot", project, sep=".")) |>
			select(feature_id, uniprot, assay)

## check ids match between pheno and prot
stopifnot(identical(pheno$patient.id, colnames(prot)))

prot_long <- prot |>
  as.data.frame() |>
  tibble::rownames_to_column("feature_id") |>
  pivot_longer(-feature_id, names_to = "patient.id", values_to = "npx") |>
  left_join(annot, by = c("feature_id"))

## ---- dist -------------------------------------------------------------
dist <- 
    prot_long |>
        ggplot(aes(x = "", y = npx)) +
		geom_violin(fill = "#F8766D", alpha = 0.5, colour = NA) +
        geom_boxplot(width = 0.15, outlier.size = 0.5, fill = "white", colour = "grey30") +
        facet_wrap(~ assay, scales = "free_y", ncol = 3) +
        labs(x = NULL, y = "NPX") +
        theme_bw(base_size = 11) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())