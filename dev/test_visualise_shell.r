## Testing the visualisation of a shell cluster in 3D
# Checking that particles are positioned and oriented as expected

library(cd)
library(rgl)
R = 8;
rho_dye = 0.1; #/nm^2
N = dye_coverage( rho_dye, R )

a = 0.5; b=a; c=1.2;

cl = cluster_shell(N, R, c, a, b, c, 'radial', 'random');
cl = cluster_shell(N, R, c, a, b, c, 'flat', 'hc');
cl = cluster_shell(N, R, c, a, b, c, 'radial', 'fibonacci');

visualise_cluster(cl)
rgl::rgl.spheres(0,0,0, cl$R0)
# cd:::sample_random(3)

cl$colours <- col2rgb(scales::hue_pal()(ncol(cl$sizes)))/255
povray_cluster(cl, "shell.pov")
