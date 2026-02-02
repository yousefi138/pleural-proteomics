# packages
packages <- c("eval.save", "ewaff", "purrr", "dplyr")
#, "aries"
lapply(packages, require, character.only=T)

# dirs
# set dirs  
dir <- paths
eval.save.dir(dir$cache)

# src the protein summary function
source(file.path(dir$scripts, "R/protein.summary.r"))

# output specs
out <- list()
#out$file.prefix <- "ewas-alspac-"

## ----access.pheno -------------------------------------------------------------
## retrieve batch info for project
batch <- eval.ret(paste("batch", project, sep="."))
batch$plate <- as.factor(batch$plateid)

# Add batch data and restrict pheno to those with proteins passing qc
pheno <- eval.ret("pheno") |>
			left_join(batch, by = c("patient.id" = "sample_id")) |>  

           ## keep only samples passing qc
            filter(excluded == FALSE)

## ----get protins -------------------------------------------------------------
prot.mat <- eval.ret(paste("prot.mat", project, sep="."))
prot <- prot.mat[,match(pheno$patient.id, colnames(prot.mat))]

## check ids match between pheno and prot
identical(pheno$patient.id, colnames(prot))

## check missingness per protein
apply(prot, 1, function(i) sum(is.na(i)))

## ----source-models -------------------------------------------------------------
source(paste0("models-", project,".r"), echo=T, max.deparse.length = 500)
models

## ----run -------------------------------------------------------------
inputs	<-
		list(model = models,
			vars = map(model.vars, ~ .x[1]),
			model.names = names(models))

eval.save({

	ret <- inputs |>
			pmap(~ {
				ret <- ewaff.sites(..1, 
						variable.of.interest = ..2, 
						methylation = prot, 
						data = pheno, 
						family="gaussian",
						method="glm",
						generate.confounders=NULL)

				sum.ret <- protein.summary(ret, 
						molecules = prot)

				output.file = file.path(dir$scripts, "docs", paste0(..3, ".html"))

				protein.report(sum.ret,
					output.file = output.file,
					author = "Paul Yousefi",
					study = paste0("Pleural proteomics analysis of ", ..3, "variable"))

				file.copy(output.file,
					file.path(dir$output, "."), overwrite = T)

				list(ret = ret,
					sum.ret = sum.ret)
				})

}, "ret", redo=T)
ret <- eval.ret("ret")


