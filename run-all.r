args <- commandArgs(trailingOnly=TRUE)

config.name <- "local"
if (length(args) > 0)
    config.name <- args[1]

paths <- config::get(config=config.name)
print(paths)

paths$scripts <- file.path(paths$project, "scripts")
paths$data <- file.path(paths$project, "data")
paths$output <- file.path(paths$project, "results")
paths$cache <- file.path(paths$project, "results", "analysis-cache")
if(!dir.exists(paths$cache)) dir.create(paths$cache)

print(paths)

## clean raw olink protein data
## in: "GB390725-RB_pleural fluid_NPX.csv"
## out: "metaboprep_export/"
##      "project_metaboprep_qc_report.html"
##      "project_metaboprep_qc_report.log"
source("proteins.r", echo=T, max.deparse.length = 500)


## clean raw pheno data
## in: "Proteomics Infection and Controls 10.11.25.xlsx"
##      "metaboprep_export/qc/data.tsv" 
## out: pheno.rda in analysis-cache
source("pheno.r", echo=T, max.deparse.length = 500)



## run analysis looking at relationship between
## methylation predicted proteins and 
## tumor vs. normal tissue type. 
## render an html summary
packages <- c("rmarkdown", "knitr")
lapply(packages, require, character.only=T)

render("analysis.rmd", 
	output_format = "all",
    output_dir = "docs")
