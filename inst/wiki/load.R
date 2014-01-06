
## ----load,echo=FALSE,warning=FALSE--------------------------------------------
library(cda)
library(rgl)
library(ggplot2)
library(reshape2)
library(plyr)
library(knitr)
#require(microbenchmark)
#require(xtable)

## ----setup,echo=FALSE----------------------------------------------------
knit_hooks$set(rgl = function(before, options, envir) {
  # if a device was opened before this chunk, close it
  if (before && rgl.cur() > 0) rgl.close()
  hook_rgl(before, options, envir)
})
rgl_annotate = function(){
  axes3d( labels = FALSE, tick = FALSE, edges=c("x", "y", "z") )
axis3d(labels = FALSE, tick = FALSE, 'x',pos=c(NA, 0, 0))
axis3d(labels = FALSE, tick = FALSE, 'y',pos=c(0, NA, 0))
axis3d(labels = FALSE, tick = FALSE, 'z',pos=c(0, 0, NA))
title3d('','','x axis','y axis','z axis')
}
theme_set(theme_minimal())

display_benchmark <- function(x, unit = "t"){
  x$time <- microbenchmark:::convert_to_unit(x$time, unit)
  res <- aggregate(time ~ expr, x, function(.x) c(mean(.x), median(.x), min(.x), max(.x)))
  res <- cbind(res$expr, as.data.frame(res$time))
  colnames(res) <- c("expr",  "mean", "median", "min", "max")
  print(xtable(res), type = 'html', html.table.attributes = '')
}