## ----globals -------------------------------------------------------------
packages <- c("dplyr") 
lapply(packages, require, character.only=T)

# set dirs  
dir <- paths
eval.save.dir(dir$cache)

## directory for metaboprep export
dir$project <-file.path(dir$output, project)
if(!dir.exists(dir$project)){
	dir.create(dir$project)
}
## ----prep.elisa -------------------------------------------------------------
## load raw clincal phenotype info 
elisa <-  data.table::fread(file.path(dir$data, 
                "il6_clean_los.txt")) 
colnames(elisa) <- colnames(elisa) |>
                make.names()|>
                tolower()
elisa <- as_tibble(elisa)				
str(elisa)

# just keep relevant data
# mostly restricting to drop out DOB and any other potentially sensitive info
out <- elisa |>
		dplyr::select(patient_id, age, sex, serum_il6, pleural_il6, 
			neutrophils, crp) |>
		rename(elisa_serum_il6 = serum_il6,
			elisa_pleural_il6 = pleural_il6)

data.table::fwrite(out, file.path(dir$project, 
                "pleural-serum-il6-elisa.csv")) 

## ----copy.olink -------------------------------------------------------------
olink <- list()
olink$pleural <- "GB390725-RB_pleural fluid_NPX.csv"
olink$plasma <- "GB390725-RB_plasma_NPX.csv"

sapply(olink, function(i) {
	file.copy(
		file.path(dir$data, i),
		file.path(dir$project, i))
})

## ----make-readme -------------------------------------------------------------
readme_content <- "# University of Bristol Pleural Proteomics README

## Contacts
- This data was prepared by Paul Yousefi (paul.yousefi@bristol.ac.uk) and transferred on 16 Feb 2026
- David Arnold (dt.arnold@bristol.ac.uk) is study lead and primary contact

## Description
This project generated 3x plates of Olink Target 96 Inflammation pannel 
data on pleural fluid samples. We also generated 1x plate of the same
assay on plasma samples from a subset of participants with pleural 
pannel data. 

Again, for a subset of these individuals, pleural and serum il6 measurements
by ELISA assay were available from a previous project.

Full protein abundance data for all 3 of these project

## Files
- GB390725-RB_pleural fluid_NPX.csv # complete pleural olink data
- GB390725-RB_plasma_NPX.csv # complete plasma olink 
- pleural-serum-il6-elisa.csv # complete elisa pleural and serum il6 data
"

writeLines(readme_content, file.path(dir$project, "README.md"))

## ----zip -------------------------------------------------------------
# Create a gzipped tar archive of the directory
tar(
  tarfile = file.path(dir$output, paste0(project, ".tar.gz")),
  files = dir$project,
  compression = "gzip"
)

