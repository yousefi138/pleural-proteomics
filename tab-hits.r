## ----fig-assoc -------------------------------------------------------------
packages <- c("eval.save","dplyr", "purrr") 
lapply(packages, require, character.only=T)

dir <- paths
eval.save.dir(dir$cache)

## ----load.data -------------------------------------------------------------
projects <- c("pleural", "plasma")
annot <- eval.ret(paste("annot", projects[1], sep="."))

## res
model <- "infect.bi.new.fulladj"

hits <- do.call(rbind,
        sapply(projects, function(project){
                ret <- eval.ret(paste("ret", project, sep="."))
                idx <- grep(model, names(ret))
                df <- ret[[idx]]$ret$table
                sig <- ret[[idx]]$sum.ret$significant.sites
                df$model_name <- project
                df$feature_id <-rownames(df)
                stopifnot(identical(rownames(df), annot$feature_id))
                df <- df |>
                      mutate(uniprot = annot$uniprot,
                          gene = annot$assay,
                          model_name = factor(model_name, 
                              levels = c("pleural", "plasma"))) 
				df <- subset(df, feature_id %in% sig) |>
					dplyr::arrange(p.value) |>
					dplyr::select(gene, n, estimate, se, p.value, p.adjust, model_name) |>
					dplyr::mutate(across(c(p.value, p.adjust),
										~ formatC(., format = "e", digits = 1))) 
        df
        }, simplify=F)
)
