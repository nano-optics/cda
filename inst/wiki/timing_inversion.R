
## @knitr load
require(cda)
require(ggplot2)
require(microbenchmark)
require(xtable)

theme_set(theme_minimal())
display_benchmark <- function(x, unit = "t"){
    x$time <- microbenchmark:::convert_to_unit(x$time, unit)
    res <- aggregate(time ~ expr, x, function(.x) c(mean(.x), median(.x), min(.x), max(.x)))
    res <- cbind(res$expr, as.data.frame(res$time))
    colnames(res) <- c("expr",  "mean", "median", "min", "max")
    print(xtable(res), type = 'html', html.table.attributes = '')
}


## @knitr cluster
wvl <- c(400, 900)
gold <- epsAu(wvl)

cl1 <- cluster_chain(N=2, pitch=500)
cl2 <- cluster_chain(N=30, pitch=500)

Angles1 <- c(0, pi/2)
Axes1 <- 'z'

Angles2 <- rep(seq(0, pi/2, length=300), 3)
Axes2 <- rep(c('x','y','z'), each=300)


## @knitr unnamed-chunk-1
comparison1 <- microbenchmark(direct = dispersion_spectrum(cl1, Angles1, Axes1, gold, invert=FALSE),
                              inversion = dispersion_spectrum(cl1, Angles1, Axes1, gold, invert=TRUE))
display_benchmark(comparison1)


## @knitr unnamed-chunk-2
comparison2 <- microbenchmark(direct = dispersion_spectrum(cl1, Angles2, Axes2, gold, invert=FALSE),
                              inversion = dispersion_spectrum(cl1, Angles2, Axes2, gold, invert=TRUE))
display_benchmark(comparison2)


## @knitr unnamed-chunk-3
comparison3 <- microbenchmark(direct = dispersion_spectrum(cl2, Angles1, Axes1, gold, invert=FALSE),
                              inversion = dispersion_spectrum(cl2, Angles1, Axes1, gold, invert=TRUE))
display_benchmark(comparison3)


