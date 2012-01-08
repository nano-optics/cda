
background { color rgb <1.0,1.0,1.0> }
global_settings {
ambient_light <0.3,0.3,0.3>
}

#include "colors.inc"
#include "textures.inc"
#include "glass.inc" 

#include "axes.pov"

light_source{<20,30,-20> color rgbf <1, 1, 1, 0.9> 
//shadowless
 fade_distance 20
    fade_power 0.5
}

#declare R = 1.0;
#declare dist = 0.8;


camera {location <dist , dist ,dist>
right <0, 4/3, 0>
 up    <0,0,1>
 look_at  <0.0 , 0.0 , 0.0>}
       

// Ground plane 
plane { y, -30*R
texture{
                   // pigment {color rgbf <0.85, 0.9, 1, 0.9> }
                    pigment {color rgbf <0.8, 0.8, 0.8, 0.9> }
                    finish{
                        reflection {1,0} 
                        ambient 1 diffuse 1  
                        specular 0.6 roughness 1/10
                    }                                        
                }     
finish {  F_Glass10 } interior {I_Glass1 } 
}
#declare Chain =	sphere {
	
	<0,0,0>, 0.002
	 no_shadow
       no_reflection
       pigment { agate }
     finish {
        ambient .1
        diffuse .1
        specular 0.9
       phong 1
        roughness .01
        metallic
        reflection {
          .8
          metallic
        }    
}
}
       // Particle
#declare Particle =	sphere {
	
	<0,0,0>, R
	 no_shadow
       //no_reflection
	pigment { BrightGold }
	normal {dents 0.0 scale 0.0}
     finish {
        ambient .1
        diffuse .01
        specular 0.8
       phong 1
        roughness .01
        metallic
        reflection {
          .8
          metallic
        }    
}
}
// Euler Angles in Bunge convention

#include "euler.pov"

//
#declare phi = 30;
#declare theta = 0;
#declare psi = -30;

//object{ Particle scale <0.05,0.03,0.03> rotate<0, 0, 0> translate <0.1, 0, 0> }
//object{ Particle scale <0.1,0.01,0.01> EulerRotate(<phi, theta , psi>) translate <0, 0, 0> }

#include "pos-helix-ell.pov"
    

//#include "positions2.pov"
object{ AxisXYZ( 0.5, 0.5, 0.5, Texture_A_Dark, Texture_A_Light)} 