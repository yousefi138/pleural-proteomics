## ----globals -------------------------------------------------------------
packages <- c("eval.save", "knitr", "corrplot", "dplyr", "ggplot2")
# "tableone", "tidyverse") 
lapply(packages, require, character.only=T)

dir <- paths
eval.save.dir(dir$cache)

## ----load.data -------------------------------------------------------------

## proteins
plasma.all <- eval.ret(paste("prot.mat", "plasma", sep="."))
pleural.all <- eval.ret(paste("prot.mat", "pleural", sep="."))

## annot
annot <- eval.ret(paste("annot", "pleural", sep=".")) |>
			select(feature_id, assay) |>
			mutate(protein = make.names(tolower(assay)))

## ----overlap -------------------------------------------------------------
## overlap
kable(addmargins(table(colnames(pleural.all) %in% colnames(plasma.all))))

## ----restrict -------------------------------------------------------------
ids <- intersect(colnames(pleural.all), colnames(plasma.all))
plasma <- plasma.all[,ids]
pleural <- pleural.all[,ids]
identical(colnames(pleural), colnames(plasma))
identical(rownames(pleural), rownames(plasma))

## ----correlation ----------------------------------------------------
## compute correlation between proteins across both tissues
pearson <- cor(t(plasma), t(pleural), 
			method = "pearson",
			use="pairwise.complete.obs")

identical(annot$feature_id, rownames(pearson))
colnames(pearson) <- rownames(pearson) <- annot$protein

## ----cor.mat -------------------------------------------------------------
corrplot(pearson, method="color", type="full", 
		xlab = "plasma proteins", ylab = "pleural proteins",
		diag=TRUE, tl.cex=0.7, title="Protein Correlations")

## ----hist -------------------------------------------------------------
cors <-  data.frame(pearson = diag(pearson))
ggplot(cors, aes(x = pearson)) +
	geom_histogram(aes(y = after_stat(density))) +
	geom_density() +
    scale_x_continuous(limits = c(-1, 1)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red")

## il6 correlation
il6 <- which(annot$protein %in% "il6")
cors[il6, ]
