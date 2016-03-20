module spherical
  implicit none

CONTAINS

  integer function cart2sphere(x, y, z, r, t, p)
    implicit none

    real*8 x, y, z 
    !	cartesian coordinates of point. Input arguments
    
    real*8 r, t, p
    !	sperical coordinates of point. Output arguments
    !	radius, theta, phi
    !	theta is [0, pi]
    !	phi is [-pi, pi]
    
    r = sqrt(x*x+y*y+z*z)
    
    If 	(r .eq. 0d0) then
       t= 0d0
       p= 0d0		
    else
       t = asin(z / r)
       
       if( (y .eq. 0d0).and.(x .eq. 0d0) ) then
          p=0
       else 
          p = atan2(y,x)
       end if
    end if
    
    cart2sphere = 1d0
  end function Cart2sphere


  real*8 function sphere2cart(x, y, z, r, t, p)
    implicit none
    
    real*8 r, t, p
    !	sperical coordinates of point. Input arguments
    !	radius, theta, phi
    !	theta is [-pi/2, pi/2]
    !	phi is [-pi, pi]
    
    real*8 x, y, z 
    !	cartesian coordinates of point. Output arguments
    
    x = r * cos(t) * cos(p)
    y = r * cos(t) * sin(p)
    z = r * sin(t)
    
    sphere2cart = 1
    
  end function Sphere2cart

end module spherical