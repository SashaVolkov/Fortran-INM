module geometry

	implicit none

CONTAINS


	subroutine distance_sphere(radius, angl1, angl2, dist)
		real(8), intent(in) :: angl1(1:2),angl2(1:2), radius
		real(8), intent(out) :: dist
		real(8) r1(1:3), r2(1:3)

		call sphere2car(radius, angl1, r1)
		call sphere2car(radius, angl2, r2)
		dist = dist_r(r1, r2)

	end subroutine


	real(8) function dist_r(r1,r2)
		real(8), dimension(1:3) :: r1, r2
		dist_r = dsqrt(sum((r2-r1)*(r2-r1)))
	end function dist_r

	subroutine sphere2car(radius, lat_lon, r_vec)

		real(8), intent(in) :: radius, lat_lon(1:2)
		real(8), intent(out) :: r_vec(1:3)

		r_vec(1) = radius * dcos(lat_lon(2)) * dcos(lat_lon(1))
		r_vec(2) = radius * dcos(lat_lon(2)) * dsin(lat_lon(1))
		r_vec(3) = radius * dsin(lat_lon(2))

	end subroutine Sphere2car


end module geometry
