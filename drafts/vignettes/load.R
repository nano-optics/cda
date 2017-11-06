
## ----load,echo=FALSE,warning=FALSE--------------------------------------------
library(cda)
library(rgl)
library(ggplot2)
library(reshape2)
library(plyr)
library(dielectric)
library(knitr)
#require(microbenchmark)
#require(xtable)

## ----setup,echo=FALSE----------------------------------------------------
knit_hooks$set(rgl = function(before, options, envir) {
  # if a device was opened before this chunk, close it
  if (before && rgl.cur() > 0) rgl.close()
  hook_rgl(before, options, envir)
})

theme_set(theme_minimal())

display_benchmark <- function(x, unit = "t"){
  x$time <- microbenchmark:::convert_to_unit(x$time, unit)
  res <- aggregate(time ~ expr, x, function(.x) c(mean(.x), median(.x), min(.x), max(.x)))
  res <- cbind(res$expr, as.data.frame(res$time))
  colnames(res) <- c("expr",  "mean", "median", "min", "max")
  print(xtable(res), type = 'html', html.table.attributes = '')
}