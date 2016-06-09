## Script to visualise a cluster in 3D
# a spherical shell with prolate particles of arbitrary size, arbitrary positions, and
# arbitrary orientations. Random colours are assigned for completeness



library(cd)

N = 125;
R = 15;
s = 2;

cl = cluster_shell(N, R, s, 1, 1, 1,  'flat', 'fibonacci');
cl$sizes = s*rbind(runif(N, 0.2, 1.2), runif(N, 0.2, 1.2), runif(N, 0.2, 1.2));
colours = scales::hue_pal()(N)

library(rgl)

visualise(cl, col=colours)

visualise(cl, "povray", out="cluster.pov")


