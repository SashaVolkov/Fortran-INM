module spherical
	implicit none
		!	radiusius, longitude, latitude
		!	longitude is [-pi/2, pi/2]
		!	latitude is [-pi, pi]
CONTAINS

	subroutine cart2sphere(x, y, z, radius, longitude, latitude)

		real(8), intent(in) :: x, y, z 
		real(8), intent(out) :: radius, longitude, latitude

		radius = sqrt(x*x+y*y+z*z)


			longitude = dasin(z / radius)
			
			if( (y == 0d0).and.(x == 0d0) ) then
				latitude=0
			else 
				latitude = datan2(y,x)
			end if

	end subroutine Cart2sphere


	subroutine sphere2cart(x, y, z, radius, longitude, latitude)

		real(8), intent(in) :: radius, longitude, latitude
		real(8), intent(out) :: x, y, z 

		x = radius * dcos(longitude) * dcos(latitude)
		y = radius * dcos(longitude) * dsin(latitude)
		z = radius * dsin(longitude)

	end subroutine Sphere2cart

end module spherical