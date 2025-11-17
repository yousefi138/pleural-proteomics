## ----globals -------------------------------------------------------------
packages <- c("readxl", "eval.save", "dplyr") 
lapply(packages, require, character.only=T)

# set dirs  
dir <- paths
eval.save.dir(dir$cache)

## ----load.data -------------------------------------------------------------
list.files(dir$data)

## load raw clincal phenotype info 
raw <- read_excel(file.path(dir$data, 
                "Proteomics Infection and Controls 10.11.25.xlsx")) 
colnames(raw) <- colnames(raw) |>
                make.names()|>
                tolower()
str(raw)

## load the proteins for patient.ids surving qc
prot <- data.table::fread(file.path(dir$output,
                "metaboprep_export/qc/data.tsv"))$sample_id

## retrieve batch info
batch <- eval.ret("batch")

# make a factor variable for plate
batch$plate <- as.factor(batch$plateid)

## ----make.pheno -------------------------------------------------------------
pheno <- raw |>    
            left_join(batch, by = c("patient.id" = "sample_id")) |>        
            mutate(female = sign(sex == "Female"),
					age  = age.at.enrollment) |>
            mutate(infect.fct = {
                factor(
                    ifelse(grepl("SPE", raw$final.diagnosis.1), 1, 
                        ifelse(grepl("CPPE", raw$final.diagnosis.1), 2,
                            ifelse(grepl("utrue|bacterial", raw$final.diagnosis.1), 3, NA))),
					levels = c(1, 2, 3),
					labels = c("control", "inter", "case"))                        
                        }) |>
            mutate(infect.num = as.numeric(infect.fct)) |>
            mutate(infect.bi = sign(infect.num>2)) |>
			relocate(patient.id, age, female, infect.fct, infect.num, infect.bi) |>
            
            ## keep only samples passing qc
            filter(patient.id %in% prot) |>
            eval.save("pheno", redo=T)            
pheno <- eval.ret("pheno")


table(pheno$infect.fct)                
table(pheno$infect.fct, pheno$infect.num)
table(pheno$infect.bi, pheno$infect.num)

table(is.na(pheno$infect.fct))
table(is.na(pheno$infect.num))
table(is.na(pheno$infect.bi))

## check the pheno$infect categories were created correctly
sum(grepl("SPE", pheno$final.diagnosis.1))
sum(grepl("CPPE", pheno$final.diagnosis.1))
sum(grepl("utrue|bacterial", pheno$final.diagnosis.1))



