visualise_povray <- function(cl, outfile="positions.pov"){

  if(is.null(cl$colours)){
cat(paste("object{ Particle scale <", 
          apply(round(cl$sizes, 5), 2, paste, collapse=", "), "> EulerRotate(<",
          apply(round(cl$angles, 5), 2, paste, collapse=", "), ">) translate <",
          apply(round(cl$positions, 5), 2, paste, collapse=", "), "> }", sep=""), 
    sep="\n", file=outfile, append=FALSE) 
  } else {
    cat(paste("sphere {
	            <0,0,0>, R
              scale <", 
              apply(round(cl$sizes, 5), 2, paste, collapse=", "), "> 
              EulerRotate(<",
              apply(round(cl$angles, 5), 2, paste, collapse=", "), ">) 
              translate <",
              apply(round(cl$positions, 5), 2, paste, collapse=", "), " >
              pigment {rgb <", 
              apply(round(cl$colours, 5), 2, paste, collapse=", "), "> }
              //no_shadow
              no_reflection
              //normal { bumps 1 scale 0.15}
              normal { dents 0 scale 1/5}
              finish {
              ambient .01
              diffuse .01
              specular 0.9
              phong 10
              roughness .01
              metallic
              reflection {
              .9
              metallic
              }    
              
              }
  }", sep=""), sep="\n",
        file=outfile, append=FALSE) 
    
    }

}

