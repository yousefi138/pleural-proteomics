kable_my_defaults <- function(x, ...){
	require('knitr')
	require('kableExtra')

	kable(x, digits = 4) %>%
  	kable_styling( 
  		bootstrap_options = c("striped", "hover"), 
  		full_width = F,
  		position = "left",...)
}