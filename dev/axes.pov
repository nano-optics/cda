
#declare Texture_A_Dark  = texture {
                               pigment{color rgb<1,0.25,0>}
                               finish {ambient 0.15 diffuse 0.85 phong 1}
                             }
#declare Texture_A_Light = texture { 
                               pigment{color rgb<1,1,1>}
                               finish {ambient 0.15 diffuse 0.85 phong 1}
                             }
                             
                             
#declare Texture_A_Dark2  = texture {
                               pigment{color rgb<0.2,0,1>}
                               finish {ambient 0.15 diffuse 0.85 phong 1}
                             }
#declare Texture_A_Light2 = texture { 
                               pigment{color rgb<1,1,1>}
                               finish {ambient 0.15 diffuse 0.85 phong 1}
                             }


#macro Axis_( AxisLen, Dark_Texture,Light_Texture, scaling) 
 union{
    cylinder { <0,-AxisLen,0>,<0,AxisLen,0>, 0.001*1000*scaling
                   texture{checker texture{Dark_Texture } 
                               texture{Light_Texture}
                       translate<0.1,0,0.1>*1000*scaling}
             }
    cone{<0,AxisLen,0>,0.003*1000*scaling,<0,AxisLen+0.01*1000*scaling,0>,0
          texture{Dark_Texture}
         }
     } // end of union                   
#end // of macro "Axis()"
//------------------------------------------------------------------------
#macro AxisXYZ( AxisLenX, AxisLenY, AxisLenZ, scaling)
//--------------------- drawing of 3 Axes --------------------------------
union{
#if (AxisLenX != 0)
 object { Axis_(AxisLenX,  pigment{color rgb<1, 0,0>}, pigment{color rgb<1,1,1>}, scaling)   rotate< 0,0,-90>
 no_shadow
	no_reflection
}// x-Axis
#end // of #if 
#if (AxisLenY != 0)
 object { Axis_(AxisLenY,  pigment{color rgb<0,1,0>}, pigment{color rgb<1,1,1>}, scaling)  rotate< 0,0,  0>
 no_shadow
	no_reflection
}// y-Axis
#end // of #if 
#if (AxisLenZ != 0)
 object { Axis_(AxisLenZ,  pigment{color rgb<0,0,1>}, pigment{color rgb<1,1,1>}, scaling)   rotate<90,0,  0>
 no_shadow
	no_reflection
}// z-Axis
#end // of #if 
} // end of union
#end// of macro "AxisXYZ( ... )"

