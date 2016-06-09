
// Euler Angles in Bunge convention
#macro EulerRotate(A)
   #local xp = vaxis_rotate(x,z,-A.x);
   #local yp = vaxis_rotate(y,z, A.x);
   #local ys = vaxis_rotate(yp,xp,A.y);
   #local zs = vaxis_rotate(z,xp,A.y);
   #local xt = vaxis_rotate(xp,zs,A.z);
   #local yt = vaxis_rotate(ys,zs,A.z);
   
   #local cosphi = cos(A.x); 
   #local sinphi = sin(A.x);
   #local costheta = cos(A.y); 
   #local sintheta= sin(A.y); 
   #local cospsi = cos(A.z); 
   #local sinpsi = sin(A.z); 
   
      transform {
       matrix < cosphi*costheta*cospsi - sinphi*sinpsi, sinphi*costheta*cospsi + cosphi*sinpsi, -sintheta*cospsi,
                   -cosphi*costheta*sinpsi - sinphi*cospsi, -sinphi*costheta*sinpsi + cosphi*cospsi, sintheta*sinpsi,
                    cosphi*sintheta, sinphi*sintheta, costheta,
               0,0,0 >
      }
#end