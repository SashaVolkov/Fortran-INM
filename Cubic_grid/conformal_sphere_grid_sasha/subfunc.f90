module subfunc

IMPLICIT NONE

	Private
	Public :: func

	Type func
		CONTAINS
		Procedure :: cube2sphere => cube2sphere
		Procedure :: matrix_rotation_to_top => matrix_rotation_to_top
	End Type

CONTAINS


	subroutine cube2sphere(this, x, y, z, x_edge, y_edge, r_sphere, number_edge, status)
	 ! mapping between points on edge of cubed sphere on real sphere
	 ! xi_edge - coordinate on number's edge, in [-1, 1]
	 ! r_sphere - radius of inscribed sphere
	 ! x, y, z - cartesian coordinates of point

		real(8) x, y, z
		real(8) x_edge, y_edge, r_sphere
		integer number_edge, status

		real(8) rr
		Class(func) :: this

		status = 0
		rr = r_sphere / sqrt(1 + x_edge**2 + y_edge**2)

		select case (number_edge)
			 case (1)
					x = y_edge * rr
					y = x_edge * rr
					z = - rr
			 case (2)
					x = rr
					y = x_edge * rr
					z = y_edge * rr
			 case (3)
					x = - x_edge * rr
					y = rr
					z = y_edge * rr
			 case (4)
					x = - rr
					y = - x_edge * rr
					z = y_edge * rr
			 case (5)
					x = x_edge * rr
					y = - rr
					z = y_edge * rr
			 case (6)
					x = -y_edge * rr
					y = x_edge * rr
					z = rr 
			 case default
					status = 1
		end select

	end subroutine



	subroutine matrix_rotation_to_top(this,x,y,z, rot, status)
	! computes matrix rot for rotation of cube
	! nearest to (x,y,z) vertex to (0,0,R)
	! next nearest to (2sqrt2,0,3)
		
		real(8), intent(in) :: x,y,z
		real(8), dimension(1:3,1:3), intent(out) :: rot
		integer(4), intent(out) :: status
		Class(func) :: this

		rot(1:3,1) =(/  1, -2, -1 /) / sqrt(6d0)
		rot(1:3,2) =(/ -1, 0,  -1 /) / sqrt(2d0)
		rot(1:3,3) =(/  1, 1, -1 /) / sqrt(3d0)

		rot = transpose(rot)
		status = 1
	end subroutine matrix_rotation_to_top

end module