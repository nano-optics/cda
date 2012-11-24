
lf <- list.files(pattern = "rmd$")
plyr::l_ply(lf, knitr::purl)