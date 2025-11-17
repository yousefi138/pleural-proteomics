## ----globals -------------------------------------------------------------
packages <- c("eval.save", "metaboprep") 
lapply(packages, require, character.only=T)

# set dirs  
dir <- paths
eval.save.dir(dir$cache)

## ----load.data -------------------------------------------------------------
list.files(dir$data)

raw <- read_olink(file.path(dir$data, 
                "GB390725-RB_pleural fluid_NPX.csv")) 

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
    format = "html")

## ----export -------------------------------------------------------------
export(mydata, dir$output)

# take the date out of the output dir name so pipeline doesn't break 
# when there's a new expor
dir.default <- list.files(dir$output) |>
                stringr::str_subset("metaboprep_export") |>
                stringr::str_subset("[0-9]$")

if (dir.exists(file.path(dir$output, "metaboprep_export"))) {
  unlink(file.path(dir$output, "metaboprep_export"), recursive = TRUE)
}

file.rename(file.path(dir$output, dir.default),
            file.path(dir$output, sub("_[0-9].*",  "", basename(dir.default))))

## ----save.working.protein.matrix -------------------------------------------------------------
eval.save({
    prot <- data.table::fread(file.path(dir$output,
                    "metaboprep_export/qc/data.tsv")) |>
            tibble::column_to_rownames("sample_id") |>
            as.matrix()|>
            t()
}, "prot.mat", redo=T)
pheno <- eval.ret("prot.mat")

## ----save working annot -------------------------------------------------------------
## annot
eval.save({
    annot <- data.table::fread(file.path(dir$output,
                    "metaboprep_export/qc/features.tsv"))
}, "annot", redo=T)
annot <- eval.ret("annot")

## ----save batch info -------------------------------------------------------------
## annot
eval.save({
    batch <- data.table::fread(file.path(dir$output,
                    "metaboprep_export/qc/samples.tsv"))
	colnames(batch) <- colnames(batch) |>
					make.names()|>
					tolower()
	batch
}, "batch", redo=T)
batch <- eval.ret("batch")
