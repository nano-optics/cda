library(cda)
library(dielectric)
library(ggplot2)
library(purrr)
library(dplyr)

wvl <- seq(400, 800)

data(AuJC)
AuJC$set_span(400, 800)
gold <- AuJC$predict(n=300, all.knots=TRUE)

gold <- epsAu(wvl)

medium <- sqrt(1.8)
a <- 5
pitch <- 12
R0 <- 12
Na <- 6.022e23 
prefac <- (Na/2303 * 1e-14) # from nm^2 to cm^2
# 
N <- 4
cl <- cda::cluster_helix(N=N, a, a, a,
                         R0=R0,
                         pitch=pitch, delta=pi/2, right=F)


m0 <- spectrum_oa(cl, gold, medium=medium, Nq = 100) %>%
  mutate(scaled=value*prefac*N)

library(egg)
p1 <- ggplot(subset(m0,variable=="extinction" & type == "cross-section"),
       aes(energy, scaled)) +
  # facet_wrap(~type, ncol=1,scales='free')+
  geom_line() +
  scale_y_continuous(breaks=seq(0,5e8,by=1e8))+
  labs(y=expression(sigma[ext]*" (1/Mcm)"),
       x=expression(energy*" /"*eV))  +
  theme_article()
  # geom_blank(data=data.frame(type=c("cross-section","dichroism"),
  #                            scaled=c(5.01e8,2e4),energy=2))

p2 <- ggplot(subset(m0,variable=="extinction" & type == "dichroism"),
             aes(energy, scaled)) +
  # facet_wrap(~type, ncol=1,scales='free')+
  geom_line() +
  labs(y=expression(sigma[CD]*" (1/Mcm)"),
       x=expression(energy*" /"*eV))  +
  theme_article()

g <- ggarrange(p1,symmetrise_scale(p2,'y'), ncol=1)

# ggsave('gov.pdf',g,width=6,height=8)

res <- subset(m0,variable=="extinction" & type == "cross-section")[,c(1,2,6)]
names(res)[3] <- "ext"
res2 <- subset(m0,variable=="extinction" & type == "dichroism")[,c(6)]

# write.table(cbind(res,cd=res2), 'export.txt', row.names = F)

# system("open .")

model <- function(N=4, Nq=10, ...){
cl <- cda::cluster_helix(N=N, a, a, a, 
                             R0=R0,
                             pitch=step, delta=pi/2, right=F)


spectrum_oa(cl, gold, medium=medium, Nq = Nq, 
            quadrature = 'gl',iterative = F, progress = F) %>%
  mutate(scaled=value*prefac*N)
}

# m0 <- plyr::mdply(tibble(Nq=seq(50,150,by=20)), model)
# # m0 <- plyr::mdply(tibble(Nq=seq(10,50,by=10)), model)
# ggplot(subset(m0,variable=="extinction"),
#        aes(energy, scaled, colour=factor(Nq))) +
#   facet_wrap(~type, scales='free')+
#   geom_line()


