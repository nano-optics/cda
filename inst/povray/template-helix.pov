
// Global settings
background { color rgb <1.0,1.0,1.0> }
global_settings {
	ambient_light <0.3,0.3,0.3>
}
#include "colors.inc"
#include "textures.inc"
#include "glass.inc" 

light_source{<1000,2500,-2500> color White}

//dimensions
#declare Scale = 1; //units will be in microns now
#declare R = 1 * Scale;
#declare AR1 = 1;
#declare AR2 = 1;
#declare AR3 = 1;

#declare RA = R * AR1;
#declare RB = R * AR2;
#declare RC = R * AR3 ;
#declare Pitch = 1 * Scale;
#declare Camera = 0.05;


// Ground plane 
plane { y, -100*RC
finish {  F_Glass10 }	
no_shadow
no_reflection
}


// Camera settings
camera {location <0.2* Scale*Camera , 1* Scale* Camera ,-1* Scale* Camera >
right <4/3, 0, 0>
       look_at  <0.0 , 0, 0.0>}
       
       
#declare Chain =	sphere {
	
	<0,0,0>, 0.001
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
#declare Particle =	union{ sphere {
	<0,0,0>, R
	scale <RA,RB,RC>
	pigment { BrightGold }
//no_shadow
no_reflection
	//normal { bumps 1 scale 0.15}
normal { dents 0 scale 1/5}
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
}}

#include "euler.pov"
#include "pos-helix.pov"
