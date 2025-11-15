# packages
packages <- c("eval.save", "ewaff", "purrr") 
# "tidyverse", "aries"
lapply(packages, require, character.only=T)

# dirs
# set dirs  
dir <- paths
eval.save.dir(dir$cache)

# src the protein summary function
source(file.path(dir$scripts, "pleural-proteomics/R/prot.summary.r"))

# output specs
out <- list()
#out$file.prefix <- "ewas-alspac-"

## ----access.pheno -------------------------------------------------------------
pheno <- eval.ret("pheno")

## ----get dnam -------------------------------------------------------------
# pull the aries dnam from the time points and observations
# specified in 'alspac'

# match gives the index of where in the second argument the first argument is,
# hence if you index the second argument on it's result, it will be in the same order as
# the first

eval.save({
    prot <- data.table::fread(file.path(dir$output,
                    "metaboprep_export/qc/data.tsv")) |>
            tibble::column_to_rownames("sample_id") |>
            as.matrix()|>
            t()
    prot <- prot[,match(pheno$patient.id, colnames(prot))]
}, "prot", redo=F)
prot <- eval.ret("prot")

## check ids match between pheno and prot
identical(pheno$patient.id, colnames(prot))

## check missingness per protein
apply(prot, 1, function(i) sum(is.na(i)))

## ----define models -------------------------------------------------------------
model.vars <- list("infect", "female", "age")

models <- 
	model.vars |>
		map(~{
				reformulate(c(.x), response = "methylation")
		})
names(models) <- model.vars

## ---- -------------------------------------------------------------
inputs	<-
		list(model = models,
			vars = model.vars)

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

				list(ret = ret,
					sum.ret = sum.ret)
				})

}, "ret", redo=F)
ret <- eval.ret("ret")

