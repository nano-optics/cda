
// Euler Angles in Bunge convention
#macro EulerRotate(A)
   #local xp = vaxis_rotate(x,z,-A.x);
   #local yp = vaxis_rotate(y,z, A.x);
   #local ys = vaxis_rotate(yp,xp,A.y);
   #local zs = vaxis_rotate(z,xp,A.y);
   #local xt = vaxis_rotate(xp,zs,A.z);
   #local yt = vaxis_rotate(ys,zs,A.z);
   
   #local cosphi = vdot(xp,x); 
   #local sinphi = vdot(vcross(xp,x), z);
   #local costheta = vdot(zs,z); 
   #local sintheta= vdot(vcross(zs,z), xp); 
   #local cospsi = vdot(xt,xp); 
   #local sinpsi = vdot(vcross(xt,xp), zs); 
   
      transform {
       matrix < cospsi, sinpsi, 0,
               -sinpsi, cospsi, 0,
               0, 0, 1, 0,0,0 >
      matrix < 1, 0, 0,
               0, costheta, sintheta,
               0, -sintheta, costheta, 0,0,0 >
                   	      matrix < cosphi, sinphi, 0,
               -sinphi, cosphi, 0,
               0, 0, 1, 0,0,0 >

   }
#end