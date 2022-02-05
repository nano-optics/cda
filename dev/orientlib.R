library(cda)
library(orientlib)
library(rayvertex)

cl <- cluster_helix(12, a = 15,b=15,c=30, R0 = 80, pitch=200, delta = pi/4)


#313 to 123
euler_pyr <- function(angles){
  # angles <- rev(angles)
  cf = cos(angles[1]); sf = sin(angles[1])
  ct = cos(angles[2]); st = sin(angles[2])
  cp = cos(angles[3]); sp = sin(angles[3])
  c(atan2(cf*st, ct),
    -asin(sf*st),
    atan2(cf*sf + sf*ct*cp, cf*cp - sf*ct*sp))
  
}

convert_angles <- function(angles){
  # m <- cda:::cpp_euler_active(angles[1],angles[2],angles[3])
  m <- cda:::cpp_euler_passive(angles[1],angles[2],angles[3])
  v <- eulerzyx(rotmatrix(m))
  vv <- v@x
  # vv <- v@x[c(3,1,2)]
  vv[1] <- pi - vv[1]
  # vv[2] <- pi/2 - vv[2]
  vv
}



angles <- apply(cl$angles, 2, convert_angles)


# angles <- apply(cl$angles, 2, euler_pyr)

# 
# scene = sphere_mesh(position = 0*cl$positions[,1], radius = 1, 
#                         angle=c(angles[[1]]@x), scale = cl$sizes[,1],
#                         material = material_list(diffuse = "dodgerblue",type="phong"))

scene = sphere_mesh(position = cl$positions[,1], radius = 1, 
                    angle=c(rev(angles[,1])*180/pi),
                    # angle=c(cl$angles[,1]*180/pi), 
                    # angle=c(angles[,1]*180/pi),
                    scale = cl$sizes[,1],
                    # order_rotation = c(3,2,1),
                    order_rotation = c(1,2,3),
                    # order_rotation = c(3,1,2),
                    material=material_list(diffuse='gold', type="phong",toon_levels = 3, 
                                           toon_outline_width = 0.025,
                                           ambient='gold',ambient_intensity=0.3))

for(i in 2:ncol(cl$sizes)) {
  # col = hsv(runif(1))
  scene = add_shape(scene, sphere_mesh(position = cl$positions[,i], radius = 1, 
                                       angle=c(rev(angles[,i])*180/pi),
                                       # angle=c(cl$angles[,i]*180/pi), 
                                       # angle=c(angles[,i]*180/pi),
                                       scale = cl$sizes[,i],
                                       # order_rotation = c(3,2,1),
                                       order_rotation = c(1,2,3),
                                       # order_rotation = c(3,1,2),
                                       material=material_list(diffuse='gold', type="phong",toon_levels = 3, 
                                                              toon_outline_width = 0.025,
                                                              ambient='gold'
                                                              ,ambient_intensity=0.3)))
}

scene %>% rasterize_scene(light_info = directional_light(c(0.5,0.5,0)),
                          lookfrom = c(900, 400, 1000),camera_up = c(0,0,1),fov=15,background = 'white')

# formals(rasterize_scene)
