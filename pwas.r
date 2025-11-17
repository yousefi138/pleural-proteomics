# packages
packages <- c("eval.save", "ewaff", "purrr") 
# "tidyverse", "aries"
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
pheno <- eval.ret("pheno")

## ----get protins -------------------------------------------------------------
prot <- eval.ret("prot")

## check ids match between pheno and prot
identical(pheno$patient.id, colnames(prot))

## check missingness per protein
apply(prot, 1, function(i) sum(is.na(i)))

## ----define models -------------------------------------------------------------
model.vars <- list("infect.fct", "infect.num", "infect.bi", "female", "age")

models <- 
	model.vars |>
		map(~{
				reformulate(c(.x), response = "methylation")
		})
names(models) <- model.vars

## ----run -------------------------------------------------------------
inputs	<-
		list(model = models,
			vars = model.vars,
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


