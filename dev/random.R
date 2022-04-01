library(cda)
library(rgl)
R <- 30
N <- 2001
exclusion <- 1.2*sqrt(4/N)
k <- 30
Niter <- 100
seed <- 1234
set.seed(seed)
s <- R*hc_sample(N, exclusion, k, Niter)
# 
# rgl.spheres(0,0,0,R)
# rgl.spheres(s[,1],s[,2],s[,3],0.5, col="gold")

# saveRDS(s, file="test.rds")
# previous <- readRDS('test.rds')
# all.equal(s, previous)



# P <- 300
# N <- 2*P + 1
# ii <- seq(-N, N)
# lat <- asin(2*ii/(2*N+1))
# Phi <- (1+sqrt(5))/2
# long <- 2*pi*ii/Phi


s <- R*fibonacci_sample(N)

rgl.spheres(0,0,0,R)
rgl.spheres(s[,1],s[,2],s[,3],0.5, col="gold")






