## ----globals -------------------------------------------------------------
packages <- c("eval.save", "tableone", "tidyverse", "kableExtra", "ggrepel") 
lapply(packages, require, character.only=T)

dir <- paths
eval.save.dir(dir$cache)

## ----load.data -------------------------------------------------------------
## retrieve batch info for project
batch <- eval.ret(paste("batch", project, sep="."))

## pheno restricted to those with proteins passing qc
pheno <- eval.ret("pheno") |>
			left_join(batch, by = c("patient.id" = "sample_id")) |>  

           ## keep only samples passing qc
            filter(excluded == FALSE) |>
			filter(sampleqc == "PASS")

## ret
ret <- eval.ret(paste("ret", project, sep="."))

## annot
annot <- eval.ret(paste("annot", project, sep="."))

## proteins
prot.mat <- eval.ret(paste("prot.mat", project, sep="."))

## ----ret.annot -------------------------------------------------------------
# annotate the results in ret
ret.annot <- map(ret, ~{
    identical(rownames(.x$ret$table), annot$feature_id)
    .x$ret$table <- .x$ret$table |>
            rownames_to_column("feature_id") |>
            mutate(uniprot = annot$uniprot,
                gene = annot$assay) |>
            relocate(feature_id, uniprot, gene)
    .x
})

## ----source-models -------------------------------------------------------------
source(file.path(paths$scripts, models_file), echo=T, max.deparse.length = 500)

## ----tab -------------------------------------------------------------
cont <- "age"
cat <- c("female", "plate", 
        "infect.fct", "infect.num", "infect.bi",
        "comp.out", "infect.bi.new", "final.diagnosis.1")
tab <- CreateTableOne(data = pheno, 
						vars = c(cont, cat), 
						factorVars = cat)

					#	strata = "dna.type")				
tab_mat <- print(tab, printToggle = F, noSpaces = T, showAllLevels = T)
kable(tab_mat, format = "html") |>
	kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE)

## ----models -------------------------------------------------------------
unlist(models)

## ----hits -------------------------------------------------------------
hits <- sapply(ret, function(model){model$sum.ret$significant.sites})
hit.count <- data.frame(model = names(hits),
				n.bonf.sig = sapply(hits,length))
hit.count$n.fdr.sig <- map_dbl(ret, ~ {
        sum(p.adjust(.x$ret$table$p.value, method = "BH") < 0.05)
    })				
rownames(hit.count) <- NULL
kable(hit.count)

## ----volcano -------------------------------------------------------------
# make an index of the fct models
idx <- grep("fct", names(ret))
res <- do.call(rbind, 
    mapply(function(name, model) {
        df <- model$ret$table
        df$model_name <- name
        df
    }, names(ret)[-idx], ret[-idx], SIMPLIFY=FALSE)
)

fdr.thresh <- -log10(0.05/nrow(ret$pleural.infect.fct$ret$table))
volcano <- res |>
	ggplot(aes(x = estimate, y = -log10(p.value))) +
	geom_point() +
	geom_hline(yintercept = fdr.thresh, linetype='dashed') +
	facet_wrap(~model_name, ncol = 2)
volcano

## ----volcano.annot -------------------------------------------------------------
# make a vector of model names for viz
mods <- paste(project, c("infect.num.fulladj", "infect.bi.fulladj", 
        "comp.out.fulladj", "infect.bi.new.fulladj"), sep=".")
idx <- match(mods,names(ret.annot))
res <- do.call(rbind, 
    mapply(function(name, model) {
        df <- model$ret$table
        df$model_name <- name
        df
    }, names(ret.annot)[idx], ret.annot[idx], SIMPLIFY=FALSE)
)

volcano.annot <- res |>
  ggplot(aes(x = estimate, y = -log10(p.value))) +
  geom_point() +
  geom_hline(yintercept = fdr.thresh, linetype = "dashed") +
  geom_text_repel(
    data = ~ dplyr::filter(.x, -log10(p.value) > fdr.thresh),
    aes(label = gene),
    max.overlaps = Inf,   # don't silently drop labels
    size = 3
  ) +
  facet_wrap(~model_name, ncol = 2)
volcano.annot


## ----top ---------------------------------------------------------
ret.anot <- map(ret, ~{
    identical(rownames(.x$ret$table), annot$feature_id)
    .x$ret$table <- .x$ret$table |>
            rownames_to_column("feature_id") |>
            mutate(uniprot = annot$uniprot,
                gene = annot$assay) |>
            relocate(feature_id, uniprot, gene)
    .x
})

top <- map(ret.anot, ~{
    ids <- .x$sum.ret$practical.sites
    idx <- which(.x$ret$table$feature_id %in% ids)
    .x$ret$table[idx, ] |>
        dplyr::arrange(p.value)|>
        mutate(across(contains("p."), ~format(., scientific = TRUE))) |>
        kable(digits = 3)
})

## ----il6 -------------------------------------------------------------
il6 <- map(ret.anot, ~{
    ids <- "OID00482" # olinkid for il6
    idx <- which(.x$ret$table$feature_id %in% ids)
    .x$ret$table[idx, ] |>
        mutate(across(contains("p."), ~format(., scientific = TRUE))) |>
        kable(digits = 3)
})

## ----qq -------------------------------------------------------------
map(ret, ~ .x$sum.ret$qq.plot)

