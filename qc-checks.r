## ----globals -------------------------------------------------------------
packages <- c("eval.save", "metaboprep", "tidyverse") 
lapply(packages, require, character.only=T)

dir <- paths
eval.save.dir(dir$cache)

## ----load.data -------------------------------------------------------------

## read in raw olink results file
raw <- list(
		pleural = as.data.frame(data.table::fread(
					file.path(dir$data, "GB390725-RB_pleural fluid_NPX.csv"))),
		plasma = as.data.frame(data.table::fread(
					file.path(dir$data, "GB390725-RB_plasma_NPX.csv")))
		) |>
			map(~{
				out <- as_tibble(.x)
				colnames(out) <- colnames(out) |>
									make.names()|>
									tolower()
				out
			})

## read in qc'd proteins for ids of samples that passed qc
prot <- 
	list(
		pleural = eval.ret(paste("prot.mat", "pleural", sep=".")),
		plasma = eval.ret(paste("prot.mat", "plasma", sep="."))
	)

## ----calc -------------------------------------------------------------
## restrict raw res to samples that passed qc
passed <- raw |> 
			map2(prot, ~{
				.x |>
					filter(sampleid %in% colnames(.y))
			})

## calculate qc stats
qc <- 
	passed|>
		map(~{
			.x |>
				group_by(assay)|>
				summarise(
					n = n(),
					n.belowlod = sum(belowlod, na.rm = T),
					freq.belowlod = {
						sum(belowlod, na.rm = T)/n() 
					},
					n.miss = sum(is.na(npx)),
					freq.miss = sum(is.na(npx))/n()
				)				
		}) |>
		dplyr::bind_rows(.id="sample")

## ----bar.pleural -------------------------------------------------------------
## Pleural
qc |>
	filter(sample=="pleural")|>
	ggplot(aes(y = reorder(assay, n.belowlod), x = n.belowlod)) +
	geom_col() +
	labs(y = "Protein") 

## ----bar.plasma -------------------------------------------------------------
## Plasma
qc |>
	filter(sample=="plasma")|>
	ggplot(aes(y = reorder(assay, n.belowlod), x = n.belowlod)) +
	geom_col() +
	labs(y = "Protein") 

## ----freq.pleural -------------------------------------------------------------
## Pleural
qc |>
	filter(sample=="pleural")|>
	ggplot(aes(y = reorder(assay, freq.belowlod), x = freq.belowlod)) +
	geom_col() +
	labs(y = "Protein") 

## ----freq.plasma -------------------------------------------------------------
## Plasma
qc |>
	filter(sample=="plasma")|>
	ggplot(aes(y = reorder(assay, freq.belowlod), x = freq.belowlod)) +
	geom_col() +
	labs(y = "Protein") 


## ----miss.pleural -------------------------------------------------------------
## Pleural
qc |>
	filter(sample=="pleural")|>
	ggplot(aes(y = reorder(assay, n.miss), x = n.miss)) +
	geom_col() +
	labs(y = "Protein") 

## ----miss.plasma -------------------------------------------------------------
## Plasma
qc |>
	filter(sample=="plasma")|>
	ggplot(aes(y = reorder(assay, n.miss), x = n.miss)) +
	geom_col() +
	labs(y = "Protein") 