## ----globals -------------------------------------------------------------
packages <- c("eval.save", "corrplot", "dplyr")
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
# addmargins(table(colnames(pleural.all) %in% colnames(plasma.all)))
## FALSE  TRUE   Sum 
##   175    88   263 

## ----restrict -------------------------------------------------------------
ids <- intersect(colnames(pleural.all), colnames(plasma.all))
plasma <- plasma.all[,ids]
pleural <- pleural.all[,ids]
stopifnot(identical(colnames(pleural), colnames(plasma)))
stopifnot(identical(rownames(pleural), rownames(plasma)))

## ----correlation ----------------------------------------------------
## compute correlation between proteins across both tissues
pearson <- cor(t(plasma), t(pleural), 
			method = "pearson",
			use="pairwise.complete.obs")

identical(annot$feature_id, rownames(pearson))
colnames(pearson) <- rownames(pearson) <- annot$protein

## ----cor.mat -------------------------------------------------------------
corrplot(pearson, method="color", type="full",
         diag=TRUE, tl.cex=0.6, tl.col="black", tl.srt=90,
         cl.cex=0.8, mar=c(0,0,2,0))

## ----hist -------------------------------------------------------------
cors <-  data.frame(pearson = diag(pearson))
hist <- 
	ggplot(cors, aes(x = pearson)) +
		geom_histogram(aes(y = after_stat(density))) +
		geom_density() +
		scale_x_continuous(limits = c(-1, 1)) +
		geom_vline(xintercept = 0, linetype = "dashed", color = "red")
