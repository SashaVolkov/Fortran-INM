module geometry

	implicit none
CONTAINS


	subroutine distance_spher(radius,angl1, angl2, long_dist, lat_dist)
		real(8), intent(in) :: angl1(1:2),angl2(1:2), radius
		real(8), intent(out) ::  long_dist, lat_dist

		long_dist = (angl2(1) - angl1(1))*radius
		lat_dist = (angl2(2) - angl1(2))*radius

	end subroutine


end module geometry
