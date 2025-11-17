## ----globals -------------------------------------------------------------
packages <- c("eval.save", "tableone", "tidyverse") 
lapply(packages, require, character.only=T)

dir <- paths
eval.save.dir(dir$cache)

## ----load.data -------------------------------------------------------------
## pheno
pheno <- eval.ret("pheno")

## ret
ret <- eval.ret("ret")

## annot
annot <- eval.ret("annot")

## ----pheno -------------------------------------------------------------
str(pheno)

## ----tab -------------------------------------------------------------
cont <- "age"
cat <- c("female", "infect", "final.diagnosis.1")
tab <- CreateTableOne(data = pheno, 
						vars = c(cont, cat), 
						factorVars = cat)

					#	strata = "dna.type")				
print(tab, showAllLevels = TRUE)
summary(tab)

## ----models -------------------------------------------------------------
formulae <- map(ret, ~ gsub("^methylation", "proteins", .x$ret$formula))
formulae <- do.call(rbind, formulae) 
n <- unlist(map(ret, ~ nrow(.x$ret$design)))
formulae <- data.frame(number = 1:nrow(formulae),
                        model.name = rownames(formulae),
                        formula = formulae[,1],
                        n = n)
rownames(formulae) <- NULL

kable(formulae)

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

## ----vol.lab -------------------------------------------------------------
plot <- ret.anot[c("infect.num.fulladj", "infect.bi.fulladj")] |>
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

