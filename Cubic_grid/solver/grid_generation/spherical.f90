module spherical
	implicit none
		!	radiusius, longitude, latitude
		!	longitude is [-pi, pi]
		!	latitude is [-pi/2, pi/2]
CONTAINS

	subroutine cart2sphere(x, y, z, radius, latitude, longitude)

		real(8), intent(in) :: x, y, z 
		real(8), intent(out) :: radius, longitude, latitude

		radius = sqrt(x*x+y*y+z*z)

			latitude = dasin(z / radius)
			
			if( (y == 0d0).and.(x == 0d0) ) then
				longitude = 0
			else 
				longitude = datan2(y,x)
			end if

	end subroutine Cart2sphere


	subroutine sphere2cart(x, y, z, radius, latitude, longitude)

		real(8), intent(in) :: radius, longitude, latitude
		real(8), intent(out) :: x, y, z 

		x = radius * dcos(latitude) * dcos(longitude)
		y = radius * dcos(latitude) * dsin(longitude)
		z = radius * dsin(latitude)

	end subroutine Sphere2cart

end module spherical