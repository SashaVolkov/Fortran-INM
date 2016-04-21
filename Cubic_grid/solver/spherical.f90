module spherical
	implicit none
		!	radiusius, theta, phi
		!	theta is [-pi/2, pi/2]
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
			theta = dasin(z / radius)
			
			if( (y == 0d0).and.(x == 0d0) ) then
				phi=0
			else 
				phi = datan2(y,x)
			end if
		end if
		
	end subroutine Cart2sphere


	subroutine sphere2cart(x, y, z, radius, theta, phi)

		real(8), intent(in) :: radius, theta, phi
		real(8), intent(out) :: x, y, z 

		x = radius * dcos(theta) * dcos(phi)
		y = radius * dcos(theta) * dsin(phi)
		z = radius * dsin(theta)

	end subroutine Sphere2cart

end module spherical