## ----start, message=FALSE, echo=FALSE, results='hide'--------------------
library(rgl)
library(knitr)
knitr::knit_hooks$set(rgl = hook_rgl)
knitr::read_demo("clusters", package="cda")
knitr::read_demo("dimer_cd", package="cda")

