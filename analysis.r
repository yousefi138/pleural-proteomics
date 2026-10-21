## ----globals -------------------------------------------------------------
packages <- c("eval.save", "tableone", "tidyverse", "kableExtra") 
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

## ----qq -------------------------------------------------------------
map(ret, ~ .x$sum.ret$qq.plot)

## ----volcano ----------------------------------------------------------
ret$infect.fct$sum.ret$volcano.plot <- NULL
ret$infect.fct.plate$sum.ret$volcano.plot <-  NULL
ret$infect.fct.fulladj$sum.ret$volcano.plot <- NULL
map(ret, ~ .x$sum.ret$volcano.plot)

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

## ----vol.lab -------------------------------------------------------------
xi <- paste(project, c("infect.num.fulladj", "infect.bi.fulladj", 
        "comp.out.fulladj", "infect.bi.new.fulladj"), sep=".")
plot <- ret.anot[xi] |>
            map(~.x$ret$table)

plot |>
    map(~ {
        .x |>
            ggplot(aes(x = estimate, y = -log10(p.value))) +
            geom_point() +
            geom_hline(yintercept = -log10(0.05/length(.x$p.value)), 
                linetype='dashed') +
            geom_label(aes(label = gene), 
                data = ~ filter(., -log10(p.value) > -log10(0.05/length(.x$p.value))),
                vjust = -0.5, hjust = 0.5, size = 3)
})

