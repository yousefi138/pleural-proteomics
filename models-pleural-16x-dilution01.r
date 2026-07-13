## ----define models -------------------------------------------------------------
model.vars <- list("infect.fct", "infect.num", "infect.bi",
					 "comp.out", "infect.bi.new","female", "age")
model.vars <- c(model.vars, # crude 
				list(c("infect.fct",  "female", "age"),
					c("infect.num",  "female", "age"),
					c("infect.bi",  "female", "age"),
					c("comp.out",  "female", "age"),
					c("infect.bi.new",  "female", "age"))
				)

models <- 
	model.vars |>
		map(~{
				reformulate(c(.x), response = "methylation")
		})
names(models) <- map(model.vars, ~ {
					var <- .x[1]
					if (length(.x) ==2) var <- paste0(var, ".plate")
					if (length(.x) ==3) var <- paste0(var, ".fulladj")
					var
				})
names(models) <- paste(project, names(models), sep = ".")
## > models
## $`pleural-16x-dilution01.infect.fct`
## methylation ~ infect.fct
## <environment: 0x35f231350>
## 
## $`pleural-16x-dilution01.infect.num`
## methylation ~ infect.num
## <environment: 0x35f2398c8>
## 
## $`pleural-16x-dilution01.infect.bi`
## methylation ~ infect.bi
## <environment: 0x35f4cb638>
## 
## $`pleural-16x-dilution01.comp.out`
## methylation ~ comp.out
## <environment: 0x35f4cdbc0>
## 
## $`pleural-16x-dilution01.infect.bi.new`
## methylation ~ infect.bi.new
## <environment: 0x35f4d4250>
## 
## $`pleural-16x-dilution01.female`
## methylation ~ female
## <environment: 0x35f4dc950>
## 
## $`pleural-16x-dilution01.age`
## methylation ~ age
## <environment: 0x35f4e6100>
## 
## $`pleural-16x-dilution01.infect.fct.fulladj`
## methylation ~ infect.fct + female + age
## <environment: 0x35f4ecae8>
## 
## $`pleural-16x-dilution01.infect.num.fulladj`
## methylation ~ infect.num + female + age
## <environment: 0x35f4f93b0>
## 
## $`pleural-16x-dilution01.infect.bi.fulladj`
## methylation ~ infect.bi + female + age
## <environment: 0x35f4fdc78>
## 
## $`pleural-16x-dilution01.comp.out.fulladj`
## methylation ~ comp.out + female + age
## <environment: 0x35f504340>
## 
## $`pleural-16x-dilution01.infect.bi.new.fulladj`
## methylation ~ infect.bi.new + female + age
## <environment: 0x35f50a998>