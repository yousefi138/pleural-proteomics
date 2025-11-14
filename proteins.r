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
    output_dir = dir$output)

## ----export -------------------------------------------------------------
export(mydata, dir$output)