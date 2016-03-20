module morphism

implicit none

Public :: morp

Type morp
	CONTAINS
	Procedure :: reverse_stereo_projection => reverse_stereo_projection
End Type


CONTAINS

	subroutine reverse_stereo_projection(this, x_plane, y_plane,r_sphere, x, y, z, status)
	! reverse central projection from south pole of sphere
	! with centre in (0,0,0) on tangent plane in north pole (0,0,R)
	! return: 
	!  0 if ok
	!  1 if error, r_sphere < 0
		real(8) x_plane, y_plane, r_sphere, x, y, z
		integer(4), intent(out) :: status
		Class(morp) :: this


		real(8) r_coeff
		r_coeff = r_sphere / (2d0 + x_plane**2 +y_plane**2)

		if (r_sphere < 0d0) then
			 write(*,*) "Warning, radius of sphere is less than 0:"
			 write(*,*) "radius =", r_sphere
			 status = 1
		else
			 x = 2*sqrt(2d0) * x_plane * r_coeff
			 y = 2*sqrt(2d0) * y_plane * r_coeff
			 z = (2-x_plane**2-y_plane**2) * r_coeff

			 status = 0
		end if
	end subroutine


	subroutine stereo_projection(x_plane, y_plane,r_sphere, x, y, z)
	! central projection from south pole of sphere
	! with centre in (0,0,0) on tangent plane in north pole (0,0,R)
		real(8), intent(out) :: x_plane, y_plane, r_sphere
		real(8), intent(in) :: x, y, z

		r_sphere = sqrt(x**2+y**2+z**2)
		x_plane = sqrt(2d0)*x/(r_sphere+z)
		y_plane = sqrt(2d0)*y/(r_sphere+z)
	 
	end subroutine


	integer function matrix_rotation_to_top(x,y,z, rot)
	! computes matrix rot for rotation of cube
	! nearest to (x,y,z) vertex to (0,0,R)
	! next nearest to (2sqrt2,0,3)
		
		real(8), intent(in) :: x,y,z
		real(8), dimension(1:3,1:3), intent(out) :: rot

		rot(1:3,1) =(/  1, -2, -1 /) / sqrt(6d0)
		rot(1:3,2) =(/ -1, 0,  -1 /) / sqrt(2d0)
		rot(1:3,3) =(/  1, 1, -1 /) / sqrt(3d0)

		rot = transpose(rot)
		matrix_rotation_to_top = 1
	end function matrix_rotation_to_top

end module morphism