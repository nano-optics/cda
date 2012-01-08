
library(cda)
library(ggplot2)

## prescale CD to milidegrees
cs2mol <- (6.022e23 /2303 * 1e-8) # from micron^2 to cm^2
mol2milideg <- function(c=1.7e-9, l=0.1) 1e3*32.98 * c * l


load("mean.rda")
load("single.rda")
load("polydisperse.rda")
polydisperse <- melt(res, id=c("wavelength", "type"))
polydisperse$value[polydisperse$type=="cd"] <- polydisperse$value[polydisperse$type=="cd"] * mol2milideg()
str(polydisperse)

load("integration-monodisperse.rda")
monodisperse <- melt(res[[2]], id=c("wavelength", "type"))
monodisperse$value[monodisperse$type=="cd"] <- monodisperse$value[monodisperse$type=="cd"] * mol2milideg()

## reorder data to have CD on top

polydisperse$type <- relevel(polydisperse$type, "cd")
monodisperse$type <- relevel(monodisperse$type, "cd")

## scaling factor from absorbance
fac <- max(mean.abs$value)
fac2 <- max(single$value)


p <- 
ggplot(polydisperse) + facet_grid(type~., scales="free") +
  geom_path(aes(wavelength*1e3, value/max(value)*fac*cs2mol),
            size=1.5, colour="grey50") +
  geom_hline(aes(y=y), data=data.frame(y=0, type=c("cd", "abs")))+
  scale_y_continuous("Extinction (1/Mcm)                            CD (milidegrees)", expand=c(0, 0))+
  scale_x_continuous("Wavelength (nm)") +
  geom_blank(aes(x, y), data=data.frame(x=700, y=c(-58, 29, 0, 1.19e10),
                          type=rep(c("cd", "abs"), each=2)))+
  geom_path(data=mean.abs, aes(wavelength*1e3, value*cs2mol),
            colour="black", linetype=2) +
  geom_path(data=monodisperse, aes(wavelength*1e3, value/max(value)*fac2*cs2mol),
            size=1.5, colour=2, alpha=0.2) +
  geom_path(data=single, aes(wavelength*1e3, value*cs2mol),
            colour="black", linetype=3) +
  theme_minimal()+
  scale_colour_brewer(palette="Pastel1") +
  opts(legend.position="none",
       strip.text.y=theme_blank(),
       panel.grid.major = theme_line(colour = "grey90", size = 0.2), 
        panel.grid.minor = theme_line(colour = "grey95", size = 0.5)
       )


ggsave("fig3.pdf", p)
