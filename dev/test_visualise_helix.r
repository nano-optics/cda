## Script to visualise a cluster in 3D
# a spherical shell with prolate particles of arbitrary size, arbitrary positions, and
# arbitrary orientations. Random colours are assigned for completeness



library(cd)

N = 15;
cl = cluster_helix(N, 1, 1, 1.5, pitch = 10, R0=7, delta = 3*pi/N);
#cl$sizes = s*rbind(runif(N, 0.2, 1.2), runif(N, 0.2, 1.2), runif(N, 0.2, 1.2));
colours = scales::hue_pal()(N)

library(rgl)

visualise_cluster(cl, col=colours)

# povray_cluster(cl, "cluster.pov")


