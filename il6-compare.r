## ----globals -------------------------------------------------------------
packages <- c("eval.save", "dplyr", "ggplot2", "GGally") 
lapply(packages, require, character.only=T)

# set dirs  
dir <- paths
eval.save.dir(dir$cache)

## ----load.data -------------------------------------------------------------
list.files(dir$data)

## load raw clincal phenotype info 
dat <-  data.table::fread(file.path(dir$data, 
                "il6_clean_los.txt")) 
colnames(dat) <- colnames(dat) |>
                make.names()|>
                tolower()
str(dat)

