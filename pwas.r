# packages
packages <- c("tidyverse", "eval.save", "ewaff", "aries")
lapply(packages, require, character.only=T)

# dirs
dir <- paths[c("alspac.dnam", "cache", "project")]
dir$res <- file.path(dir$project, "results")
dir$scripts <- file.path(dir$project, "scripts/dnam-anti-depressants")
eval.save.dir(dir$cache)

# output specs
out <- list()
out$file.prefix <- "ewas-alspac-"

## ----access.pheno -------------------------------------------------------------
alspac <- eval.ret("alspac")

## ----get dnam -------------------------------------------------------------
# pull the aries dnam from the time points and observations
# specified in 'alspac'

# match gives the index of where in the second argument the first argument is,
# hence if you index the second argument on it's result, it will be in the same order as
# the first

eval.save({
	meth <-  
		alspac %>% 
			imap(~ 
				aries.select(dir$alspac.dnam, 
						time.point=.y) |>
					aries.methylation() %>%

					# keep only selected samples
					.[,match(.x$Sample_Name, colnames(.))])
	
	# windsize - had to use slightly strange syntax below to get map &
	#	ewaff naming to play nice 					
	meth <-                
		meth %>%
			map(function(i){
				ewaff.handle.outliers(i, method="winsorize")[[1]]
			})                
}, "meth.alspac", redo=F)
meth <- eval.ret("meth.alspac")

## ----get annotation -------------------------------------------------------------
alspac %>% 
	map(~ table(.x$beadchip))

# get annotation to use in report
annot <- list()
annot$`15up` <- meffil::meffil.get.features(featureset = "450k") 
annot$F24 <- meffil::meffil.get.features(featureset = "epic")

# coerce annot into same cpg order as meth
stopifnot(identical(names(meth), names(annot)))
annot <- map2(annot, meth, 
			~ .x[match(rownames(.y), .x$name),]%>%
					dplyr::select(name, chromosome, position))
map2(annot, meth, ~ identical(.x$name, rownames(.y)))

# check 
all(unlist(map2(annot, meth, ~ identical(.x$name, rownames(.y)))))

## ----define models -------------------------------------------------------------
model.vars <- 
			list(`15up.yph7305` = 
				list(
					x="yph7305",
					xs = c("age", "female",
							"Bcell", "CD4T", 
							"CD8T", "Eos", 
							"Mono", "Neu", "NK"),
					y = "methylation",
					name = "15up")
			)
model.vars$F24.yph7305 <- model.vars$`15up.yph7305`
model.vars$F24.yph7305$name <- "F24"


models <- 
	model.vars %>%
		map(~{
				reformulate(c(.x$x, .x$xs), response = .x$y)
		})

## ---- -------------------------------------------------------------
#check that the alspac, methylation, and model specs all have the
# same time_point ordering
identical(names(meth), names(alspac))
all(sapply(model.vars, function(i) pluck(i$name)) == names(meth))

# Use Saffari et al. 2017 sig.threshold - default is bonf
parameters = ewaff.report.parameters(sig.threshold = 3.6e-8)

inputs	<-
		list(model = models,
			vars = model.vars,
			names = names(meth),
			model.names = names(models))

eval.save({
	ret <- inputs %>% 
				pmap(~ {
					set.seed(837)	
						ret <- ewaff.sites(..1, 
								variable.of.interest = ..2$x, 
								methylation = meth[[..3]], 
								data = alspac[[..3]], 
								family="gaussian",
								#method="limma",
								#generate.confounders=NULL)
								method="glm",
								generate.confounders="sva",
								n.confounders = 10)
						
						sum.ret <- ewaff.summary(ret, 
								chr = annot[[..3]]$chromosome, 
								pos = annot[[..3]]$position, 
								methylation = meth[[..3]],
								parameters = parameters)
					
						ewaff.report(sum.ret,
							output.file = file.path(dir$res, paste0(out$file.prefix, ..4, ".html")),
							author = "Paul Yousefi",
							study = paste0(out$file.prefix, ..4))

						file.copy(file.path(dir$res, paste0(out$file.prefix, ..4, ".html")),
							file.path(dir$scripts, "reports", "."), overwrite = T)
						
						file.copy(file.path(dir$res, paste0(out$file.prefix, ..4, ".md")),
							file.path(dir$scripts, "reports", "."), overwrite = T)

					list(ret = ret,
						sum.ret = sum.ret)
				})
}, "ret.alspac", redo=F)
ret <- eval.ret("ret.alspac")

## ----export -------------------------------------------------------------
ret %>%
	imap(~ write_rds(.x, file.path(dir$res, paste0(out$file.prefix, .y, ".rds"))))

ret %>%
	imap(~ .x$ret$table %>%
				data.table::fwrite(file.path(dir$res, paste0(out$file.prefix, .y, ".csv"))))

## ----cleanup -------------------------------------------------------------
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
