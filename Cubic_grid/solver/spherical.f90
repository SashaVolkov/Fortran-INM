module spherical
	implicit none
		!	radiusius, theta, phi
		!	theta is [0, pi]
		!	phi is [-pi, pi]
CONTAINS

	subroutine cart2sphere(x, y, z, radius, theta, phi)

		real(8), intent(in) :: x, y, z 
		real(8), intent(out) :: radius, theta, phi

		radius = sqrt(x*x+y*y+z*z)

		if(radius == 0d0) then
			theta= 0d0
			phi= 0d0
		else
			theta = asin(z / radius)
			
			if( (y == 0d0).and.(x == 0d0) ) then
				phi=0
			else 
				phi = atan2(y,x)
			end if
		end if
		
	end subroutine Cart2sphere


	subroutine sphere2cart(x, y, z, radius, theta, phi)

		real(8), intent(in) :: radius, theta, phi
		real(8), intent(out) :: x, y, z 

		x = radius * cos(theta) * cos(phi)
		y = radius * cos(theta) * sin(phi)
		z = radius * sin(theta)

	end subroutine Sphere2cart

end module spherical