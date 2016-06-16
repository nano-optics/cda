Overview of the cda package
=================================================================
<!-- 
%\VignetteEngine{knitr::rmarkdown}
%\VignetteIndexEntry{overview}
-->




```r
1 + 1
```

```
## [1] 2
```

```r
lm(y ~ x, data = data.frame(x = 1:10, y = rnorm(10)))
```

```
## 
## Call:
## lm(formula = y ~ x, data = data.frame(x = 1:10, y = rnorm(10)))
## 
## Coefficients:
## (Intercept)            x  
##     -0.5374       0.1135
```


```r
# dielectric function
gold <- epsAu(seq(400, 700))
cl <- cluster_dimer(d = 50, a=30, b=10, c=10, dihedral = pi/4)
gold <- epsAu(seq(500, 800))
circular <- spectrum_oa(cl, gold, Nsca=10)
```

```
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
## 
## warning: solve(): refinement and/or equilibration not done due to crippled LAPACK
```


```r
p <- ggplot(circular, aes(wavelength, value, color=variable)) + 
  facet_grid(type~., scales="free") + geom_line() + 
  labs(x = "wavelength /nm", y = expression(sigma/nm^2), 
       colour = "variable")

p
```

![](overview_files/figure-html/plot-1.png)<!-- -->
