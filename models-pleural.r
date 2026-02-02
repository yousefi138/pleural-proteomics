## ----define models -------------------------------------------------------------
model.vars <- list("infect.fct", "infect.num", "infect.bi",
					 "comp.out", "infect.bi.new","female", "age")
model.vars <- c(model.vars, # crude 
				map(model.vars, ~c(.x, "plate")), # batch adjusted
				list(c("infect.fct",  "female", "age", "plate"),
					c("infect.num",  "female", "age", "plate"),
					c("infect.bi",  "female", "age", "plate"),
					c("comp.out",  "female", "age", "plate"),
					c("infect.bi.new",  "female", "age", "plate"))
				)

models <- 
	model.vars |>
		map(~{
				reformulate(c(.x), response = "methylation")
		})
names(models) <- map(model.vars, ~ {
					var <- .x[1]
					if (length(.x) ==2) var <- paste0(var, ".plate")
					if (length(.x) ==4) var <- paste0(var, ".fulladj")
					var
				})
names(models) <- paste(project, names(models), sep = ".")

## > models
## $pleural.infect.fct
## methylation ~ infect.fct
## <environment: 0x155e98ed8>
## 
## $pleural.infect.num
## methylation ~ infect.num
## <environment: 0x155ea3d58>
## 
## $pleural.infect.bi
## methylation ~ infect.bi
## <environment: 0x1570f22b0>
## 
## $pleural.comp.out
## methylation ~ comp.out
## <environment: 0x1570fb548>
## 
## $pleural.infect.bi.new
## methylation ~ infect.bi.new
## <environment: 0x1570fe8d0>
## 
## $pleural.female
## methylation ~ female
## <environment: 0x157105b68>
## 
## $pleural.age
## methylation ~ age
## <environment: 0x15710aef0>
## 
## $pleural.infect.fct.plate
## methylation ~ infect.fct + plate
## <environment: 0x15710e278>
## 
## $pleural.infect.num.plate
## methylation ~ infect.num + plate
## <environment: 0x1571153c0>
## 
## $pleural.infect.bi.plate
## methylation ~ infect.bi + plate
## <environment: 0x1571185f8>
## 
## $pleural.comp.out.plate
## methylation ~ comp.out + plate
## <environment: 0x157121740>
## 
## $pleural.infect.bi.new.plate
## methylation ~ infect.bi.new + plate
## <environment: 0x155f1e378>
## 
## $pleural.female.plate
## methylation ~ female + plate
## <environment: 0x155f254c0>
## 
## $pleural.age.plate
## methylation ~ age + plate
## <environment: 0x155f286f8>
## 
## $pleural.infect.fct.fulladj
## methylation ~ infect.fct + female + age + plate
## <environment: 0x155f31840>
## 
## $pleural.infect.num.fulladj
## methylation ~ infect.num + female + age + plate
## <environment: 0x155f347d8>
## 
## $pleural.infect.bi.fulladj
## methylation ~ infect.bi + female + age + plate
## <environment: 0x155f3b680>
## 
## $pleural.comp.out.fulladj
## methylation ~ comp.out + female + age + plate
## <environment: 0x155f3e618>
## 
## $pleural.infect.bi.new.fulladj
## methylation ~ infect.bi.new + female + age + plate
## <environment: 0x155f474c0>