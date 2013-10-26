
## @knitr load,message=FALSE, echo=1:4
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


## @knitr cluster, tidy=FALSE, fig.path='averaging-'
gold <- epsAu(seq(400, 900))

cl <- cluster_dimer(d=100, 
              dihedral=10*pi/180, alpha1=20*pi/180, alpha2=0,
              a=35, b=12)


cl2 <- cluster_helix(4, pitch=1000)


# # achiral cluster (plane of symmetry)
# cl2 <- cluster_dimer(d=100, 
#               dihedral=0*pi/180, alpha1=20*pi/180, alpha2=0,
#               a=35, b=12)

# achiral cluster (plane of symmetry) but in rotated frame
cl3 <- cl2
cl3$r <- t(cda$euler(pi/12,pi/5,pi/3) %*% t(cl2$r))
# rgl.ellipsoids(cl3$r, cl3$sizes, cl3$angles, col="gold")

## @knitr comparison, tidy=FALSE, fig.path='averaging-'
params <- expand.grid(N=c(10, 100, 1000), iterative=FALSE,
                       averaging=c("grid", "GL", "QMC", "cheap"),
                       stringsAsFactors=FALSE)

test <- circular_dichroism_spectrum(cluster=cl, material=gold, N=1000)
comparison <- mdply(params, circular_dichroism_spectrum, cluster=cl, material=gold)

p <- 
  ggplot(subset(comparison, type == "CD" & variable == "extinction")) + 
  facet_grid(averaging~.)+
  geom_path(aes(wavelength, value, colour=factor(N), group=N))+
  labs(y=expression(sigma[ext]/nm^2),
       x=expression(wavelength/nm), colour=expression(N))

p + geom_line(aes(wavelength, value),
              data=subset(test, type == "CD" & variable == "extinction"))

# 
## @knitr achiral, tidy=FALSE, fig.path='averaging-'
params <- expand.grid(N=c(100, 500, 800),
                       averaging=c("grid", "GL", "QMC", "cheap"),
                       stringsAsFactors=FALSE)

test <- circular_dichroism_spectrum(cluster=cl2, material=gold, precision=1e-6,
                                     N=100, averaging="QMC")
comparison <- mdply(params, circular_dichroism_spectrum, cluster=cl2, material=gold,
                    iterative=FALSE)

p <- 
  ggplot(subset(comparison, type == "CD" & variable == "extinction")) + 
  facet_grid(averaging~.)+
  geom_path(aes(wavelength, value, colour=factor(N), group=N))+
  labs(y=expression(sigma[ext]/nm^2),
       x=expression(wavelength/nm), colour=expression(N))

p + geom_line(aes(wavelength, value),
              data=subset(test, type == "CD" & variable == "extinction"))
# 
# 

## @knitr achiral, tidy=FALSE, fig.path='averaging-'
params <- expand.grid(N=c(100, 500, 800),
                      averaging=c("grid", "GL", "QMC"),
                      stringsAsFactors=FALSE)

test <- circular_dichroism_spectrum(cluster=cl2, material=gold, precision=1e-4,
                                     N=100, progress=TRUE)
comparison <- mdply(params, circular_dichroism_spectrum, cluster=cl3, material=gold)

p <- 
  ggplot(subset(comparison, type == "CD" & variable == "extinction")) + 
  facet_grid(averaging~.)+
  geom_path(aes(wavelength, value, colour=factor(N), group=N))+
  labs(y=expression(sigma[ext]/nm^2),
       x=expression(wavelength/nm), colour=expression(N))

p + geom_line(aes(wavelength, value),
              data=subset(test, type == "CD" & variable == "extinction"))
# 
# 
