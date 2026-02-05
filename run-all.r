args <- commandArgs(trailingOnly=TRUE)

config.name <- "local"
if (length(args) > 0)
    config.name <- args[1]

paths <- config::get(config=config.name)
print(paths)

paths$scripts <- file.path(paths$project, "scripts/repo","pleural-proteomics")
paths$data <- file.path(paths$project, "data")
paths$output <- file.path(paths$project, "results")
paths$cache <- file.path(paths$project, "results", "analysis-cache")
if(!dir.exists(paths$cache)) dir.create(paths$cache)
print(paths)

## clean raw olink PLEURAL protein data
## in: "GB390725-RB_pleural fluid_NPX.csv"
## out: "pleural/metaboprep_export/.",
##      "pleural_metaboprep_qc_report.html",
##      "pleural_metaboprep_qc_report.log",
##      eval.ret("prot.mat.pleural"), 
##      eval.ret("annot.pleural"), 
##      eval.ret("batch.pleural")
project <- "pleural"
file <- "GB390725-RB_pleural fluid_NPX.csv"
source("proteins.r", echo=T, max.deparse.length = 500)

## clean raw olink PLASMA protein data
## in: "GB390725-RB_plasma_NPX.csv"
## out: "plasma/metaboprep_export/"
##      "plasma_metaboprep_qc_report.html"
##      "plasma_metaboprep_qc_report.log"
##      eval.ret("prot.mat.plasma"), 
##      eval.ret("annot.plasma"), 
##      eval.ret("batch.plasma")
project <- "plasma"
file <- "GB390725-RB_plasma_NPX.csv"
source("proteins.r", echo=T, max.deparse.length = 500)

## clean raw pheno data
## in: "Proteomics Infection and Controls 10.11.25.xlsx"
## out: pheno.rda in analysis-cache i.e. eval.ret("pheno")
source("pheno.r", echo=T, max.deparse.length = 500)

## run pwas
## in: eval.ret("pheno")
##      eval.ret(paste("prot.mat", project, sep="."))
##      eval.ret(paste("batch", project, sep="."))
##      protein.summary.r, report.rmd
## out: rendered output in docs/ for each model run 
project <- "pleural"
source("pwas.r", echo=T, max.deparse.length = 500)

## run pwas
## in: eval.ret("pheno")
##      eval.ret(paste("prot.mat", project, sep="."))
##      eval.ret(paste("batch", project, sep="."))
##      protein.summary.r, report.rmd
## out: rendered output in docs/ for each model run 
project <- "plasma"
source("pwas.r", echo=T, max.deparse.length = 500)

## summarise pwas associations
##  in: eval.ret(paste("ret", project, sep="."))
## 		eval.ret("pheno")
## 		eval.ret(paste("batch", project, sep="."))
## 		annot <- eval.ret(paste("annot", project, sep="."))
##		prot.mat <- eval.ret(paste("ret", project, sep="."))
##	out: pleural-analysis.html
project <- "pleural"
output <- paste(project, "analysis.html", sep="-")
packages <- c("rmarkdown", "knitr")
lapply(packages, require, character.only=T)
source(file.path(paths$scripts, "R/kable_my_defaults.r")) # couldn't figure out how to render in the doc
render("analysis.rmd", output_file = output, output_format = "html_document", output_dir = "docs")

## summarise pwas associations
##  in: eval.ret(paste("ret", project, sep="."))
## 		eval.ret("pheno")
## 		eval.ret(paste("batch", project, sep="."))
## 		annot <- eval.ret(paste("annot", project, sep="."))
##		prot.mat <- eval.ret(paste("ret", project, sep="."))
##	out: plasma-analysis.html
project <- "plasma"
output <- paste(project, "analysis.html", sep="-")
packages <- c("rmarkdown", "knitr")
lapply(packages, require, character.only=T)
render("analysis.rmd", output_file = output, output_format = "html_document", output_dir = "docs")

## Compare pleural and plasma protein abundances
##  in: eval.ret(paste("prot.mat", "plasma", sep="."))
## 		eval.ret(paste("prot.mat", "pleural", sep="."))
## 		eval.ret(paste("annot", "pleural", sep=".")) 
##	out: tissue-compare.html
packages <- c("rmarkdown", "knitr")
lapply(packages, require, character.only=T)
render("tissue-compare.rmd", output_format = "html_document", output_dir = "docs")

## Compare pleural and plasma IL6 abundances
##  in: "il6_clean_los.txt" ELISA data
##		eval.ret(paste("prot.mat", "plasma", sep="."))
## 		eval.ret(paste("prot.mat", "pleural", sep="."))
##		eval.ret(paste("annot", "pleural", sep="."))
##	out: il6-compare.html
packages <- c("rmarkdown", "knitr")
lapply(packages, require, character.only=T)
render("il6-compare.rmd", output_format = "html_document", output_dir = "docs")
