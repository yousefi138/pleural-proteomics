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
					if (length(.x) ==3) var <- paste0(var, ".fulladj")
					var
				})
names(models) <- paste(project, names(models), sep = ".")

## > models
## $plasma.infect.fct
## methylation ~ infect.fct
## <environment: 0x156e22dd0>
## 
## $plasma.infect.num
## methylation ~ infect.num
## <environment: 0x156e280e8>
## 
## $plasma.infect.bi
## methylation ~ infect.bi
## <environment: 0x156f15ba8>
## 
## $plasma.comp.out
## methylation ~ comp.out
## <environment: 0x156f1ed98>
## 
## $plasma.infect.bi.new
## methylation ~ infect.bi.new
## <environment: 0x156f220e8>
## 
## $plasma.female
## methylation ~ female
## <environment: 0x156f25470>
## 
## $plasma.age
## methylation ~ age
## <environment: 0x156f2c6d0>
## 
## $plasma.infect.fct.fulladj
## methylation ~ infect.fct + female + age
## <environment: 0x156f31978>
## 
## $plasma.infect.num.fulladj
## methylation ~ infect.num + female + age
## <environment: 0x156f36778>
## 
## $plasma.infect.bi.fulladj
## methylation ~ infect.bi + female + age
## <environment: 0x156f39828>
## 
## $plasma.comp.out.fulladj
## methylation ~ comp.out + female + age
## <environment: 0x156f40698>
## 
## $plasma.infect.bi.new.fulladj
## methylation ~ infect.bi.new + female + age
## <environment: 0x1575fafd0>