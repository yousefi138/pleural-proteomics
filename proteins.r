## ----globals -------------------------------------------------------------
packages <- c("eval.save", "metaboprep") 
lapply(packages, require, character.only=T)

# set dirs  
dir <- paths
eval.save.dir(dir$cache)

## directory for metaboprep export
dir$project <-file.path(dir$output, project)
if(!dir.exists(dir$project)){
    dir.create(dir$project)
}

## ----load.data -------------------------------------------------------------
list.files(dir$data)

raw <- read_olink(file.path(dir$data, file)) 

mydata <- Metaboprep(data     = raw$data, 
                     features = raw$features, 
                     samples  = raw$samples)

mydata <- mydata |> quality_control( source_layer        = "input", 
                                     sample_missingness  = 0.2, 
                                     feature_missingness = 0.2, 
                                     total_peak_area_sd  = 5, 
                                     outlier_udist       = 5, 
                                     outlier_treatment   = "leave_be", 
                                     winsorize_quantile  = 1.0, 
                                     tree_cut_height     = 0.5, 
                                     pc_outlier_sd       = 5, 
                                     sample_ids          = NULL, 
                                     feature_ids         = NULL)
summary(mydata)                         

## ----report -------------------------------------------------------------
generate_report(mydata, 
    output_dir = dir$output,
    project = project,
    format = "html")

## ----export -------------------------------------------------------------
export(mydata, dir$project)

# take the date out of the output dir name so pipeline doesn't break 
# when there's a new expor
dir.default <- list.files(dir$project) |>
                stringr::str_subset("metaboprep_export") |>
                stringr::str_subset("[0-9]$")

if (dir.exists(file.path(dir$project, "metaboprep_export"))) {
  unlink(file.path(dir$project, "metaboprep_export"), recursive = TRUE)
}

file.rename(file.path(dir$project, dir.default),
            file.path(dir$project, sub("_[0-9].*",  "", basename(dir.default))))

## ----save.working.protein.matrix -------------------------------------------------------------
prot.project <- paste("prot.mat", project, sep=".")
eval.save({
    prot <- data.table::fread(file.path(dir$project,
                    "metaboprep_export/qc/data.tsv")) |>
            tibble::column_to_rownames("sample_id") |>
            as.matrix()|>
            t()
}, prot.project, redo=T)
prot <- eval.ret(prot.project)

## ----save working annot -------------------------------------------------------------
## annot
annot.project <- paste("annot", project, sep=".")
eval.save({
    annot <- data.table::fread(file.path(dir$project,
                    "metaboprep_export/qc/features.tsv"))
	colnames(annot) <- colnames(annot) |>
					make.names()|>
					tolower()
	annot                    
}, annot.project, redo=T)
annot <- eval.ret(annot.project)

## ----save batch info -------------------------------------------------------------
## batch
batch.project <- paste("batch", project, sep=".")
eval.save({
    batch <- data.table::fread(file.path(dir$project,
                    "metaboprep_export/qc/samples.tsv"))
	colnames(batch) <- colnames(batch) |>
					make.names()|>
					tolower()
	batch
}, batch.project, redo=T)
batch <- eval.ret(batch.project)
