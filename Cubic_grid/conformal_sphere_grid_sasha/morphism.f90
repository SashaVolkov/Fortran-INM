module morphism

implicit none
CONTAINS

	integer function reverse_stereo_projection(x_plane, y_plane,r_sphere, x, y, z)
	! reverse central projection from south pole of sphere
	! with centre in (0,0,0) on tangent plane in north pole (0,0,R)
	! return: 
	!  0 if ok
	!  1 if error, r_sphere < 0
		real(8) x_plane, y_plane, r_sphere, x, y, z

		real(8) r_coeff
		r_coeff = r_sphere / (2d0 + x_plane**2 +y_plane**2)

		if (r_sphere < 0d0) then
			 write(*,*) "Warning, radius of sphere is less than 0:"
			 write(*,*) "radius =", r_sphere
			 reverse_stereo_projection = 1
		else
			 x = 2*sqrt(2d0) * x_plane * r_coeff
			 y = 2*sqrt(2d0) * y_plane * r_coeff
			 z = (2-x_plane**2-y_plane**2) * r_coeff

			 reverse_stereo_projection = 0
		end if
	end function


	subroutine stereo_projection(x_plane, y_plane,r_sphere, x, y, z)
	! central projection from south pole of sphere
	! with centre in (0,0,0) on tangent plane in north pole (0,0,R)
		real(8) x_plane, y_plane, r_sphere, x, y, z

		r_sphere = sqrt(x**2+y**2+z**2)
		x_plane = sqrt(2d0)*x/(r_sphere+z)
		y_plane = sqrt(2d0)*y/(r_sphere+z)
	 
	end subroutine


	integer function matrix_rotation_to_top(x,y,z, rot)
	! computes matrix rot for rotation of cube
	! nearest to (x,y,z) vertex to (0,0,R)
	! next nearest to (2sqrt2,0,3)
		
		real(8), intent(in) :: x,y,z
		real(8), dimension(1:3,1:3), intent(inout) :: rot

		rot(1:3,1) =(/  1, -2, -1 /) / sqrt(6d0)
		rot(1:3,2) =(/ -1, 0,  -1 /) / sqrt(2d0)
		rot(1:3,3) =(/  1, 1, -1 /) / sqrt(3d0)

		rot = transpose(rot)
		matrix_rotation_to_top = 1
	end function matrix_rotation_to_top

end module morphism