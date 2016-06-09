## Script to visualise a dimer in 3D
# a dimer along z axis, 2 identical particles, dihedral angle pi/4

library(cd)
d = 50;
a = 50; b=20; c=12;
dihedral = pi/2;

cl = cluster_dimer(d, a,b,c, dihedral,0,0);
colours = rgb(rbind(c(1, 0, 0),  c(0, 1, 0)))

library(rgl)
palette(RColorBrewer::brewer.pal(7,"Set1"))


visualise_cluster(cl, col=colours)
povray_cluster(cl, "dimer.pov")



