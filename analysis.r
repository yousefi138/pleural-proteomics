## ----globals -------------------------------------------------------------
packages <- c("eval.save", "tableone", "purrr") 
lapply(packages, require, character.only=T)

dir <- paths
eval.save.dir(dir$cache)

## ----load.data -------------------------------------------------------------


## ----pheno -------------------------------------------------------------
pheno <- eval.ret("pheno")
str(pheno)

## ----tab -------------------------------------------------------------
cont <- "age"
cat <- c("female", "infect", "final.diagnosis.1")
tab <- CreateTableOne(data = pheno, 
						vars = c(cont, cat), 
						factorVars = cat)

					#	strata = "dna.type")				
print(tab, showAllLevels = TRUE)
summary(tab)

## ----models -------------------------------------------------------------
ret <- eval.ret("ret")

formulae <- map(ret, ~ gsub("^methylation", "proteins", .x$ret$formula))
formulae <- do.call(rbind, formulae) 
formulae <- data.frame(number = 1:nrow(formulae),
                        model.name = rownames(formulae),
                        formula = formulae[,1])
rownames(formulae) <- NULL

kable(formulae)

## ----qq -------------------------------------------------------------
map(ret, ~ .x$sum.ret$qq.plot)

## ----volcano ----------------------------------------------------------
ret$infect$sum.ret$volcano.plot <- NULL
map(ret, ~ .x$sum.ret$volcano.plot)