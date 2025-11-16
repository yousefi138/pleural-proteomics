## ----globals -------------------------------------------------------------
packages <- c("eval.save") 
lapply(packages, require, character.only=T)

dir <- paths
eval.save.dir(dir$cache)

## ----load.data -------------------------------------------------------------

# load pwas results
pheno <- eval.ret("pheno")
ret <- eval.ret("ret")

str(pheno)