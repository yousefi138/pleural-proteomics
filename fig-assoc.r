## ----fig-assoc -------------------------------------------------------------
packages <- c("eval.save", "patchwork","dplyr", "purrr", "ggplot2", "ggrepel") 
lapply(packages, require, character.only=T)

dir <- paths
eval.save.dir(dir$cache)

## ----load.data -------------------------------------------------------------
projects <- c("pleural", "plasma")
annot <- eval.ret(paste("annot", projects[1], sep="."))

## res
model <- "infect.bi.new.fulladj"

res <- do.call(rbind, 
        sapply(projects, function(project){
                ret <- eval.ret(paste("ret", project, sep="."))
                idx <- grep(model, names(ret))
                df <- ret[[idx]]$ret$table
                df$model_name <- project
                df$feature_id <-rownames(df)
                stopifnot(identical(rownames(df), annot$feature_id))
                df |>
                    mutate(uniprot = annot$uniprot,
                        gene = annot$assay,
                        model_name = factor(model_name, 
                            levels = c("pleural", "plasma"))) 
        }, simplify=F)
)

bonf.thresh <- -log10(0.05/length(unique(res$uniprot)))

volcanos <- res |>
  ggplot(aes(x = estimate, y = -log10(p.value))) +
  geom_point() +
  geom_hline(yintercept = bonf.thresh, linetype = "dashed") +
  geom_text_repel(
    data = ~ dplyr::filter(.x, -log10(p.value) > bonf.thresh),
    aes(label = gene),
    max.overlaps = Inf,   # don't silently drop labels
    size = 3
  ) +
  facet_wrap(~model_name, ncol = 2)

  
## ---- -------------------------------------------------------------
a <- res[res$model_name == projects[1], ]
b <- res[res$model_name == projects[2], ]

m <- merge(a[, c("uniprot", "gene", "estimate")],
           b[, c("uniprot", "estimate")],
           by = "uniprot", suffixes = paste0(".", projects))

fit <- lm(estimate.plasma ~ estimate.pleural, data = m)
s   <- summary(fit)

r2 <- s$r.squared
p  <- s$coefficients["estimate.pleural", "Pr(>|t|)"]

lab <- sprintf("R² = %.3f\np = %.2g", r2, p)

scatter <- 
	ggplot(m, aes(x = estimate.pleural, y = estimate.plasma)) +
	geom_point() +
	geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
	labs(x = paste(projects[1], "coefficient"),
		y = paste(projects[2], "coefficient")) +
	coord_equal() +
	annotate("text", x = -Inf, y = Inf, label = lab,
			hjust = -0.1, vjust = 1.2)