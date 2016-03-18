integer function cube2sphere(x, y, z, x_edge, y_edge, r_sphere, number_edge)
 ! mapping between points on edge of cubed sphere on real sphere
 ! xi_edge - coordinate on number's edge, in [-1, 1]
 ! r_sphere - radius of inscribed sphere
 ! x, y, z - cartesian coordinates of point
	implicit none
	real(8) x, y, z
	real(8) x_edge, y_edge, r_sphere
	integer number_edge

	real(8) rr

	cube2sphere = 0
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
				cube2sphere = 1
	end select

end function



integer function reverse_stereo_projection(x_plane, y_plane,r_sphere, x, y, z)
! reverse central projection from south pole of sphere
! with centre in (0,0,0) on tangent plane in north pole (0,0,R)
! return: 
!  0 if ok
!  1 if error, r_sphere < 0
	implicit none
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




integer function stereo_projection(x_plane, y_plane,r_sphere, x, y, z)
! central projection from south pole of sphere
! with centre in (0,0,0) on tangent plane in north pole (0,0,R)
	implicit none
	real(8) x_plane, y_plane, r_sphere, x, y, z

	r_sphere = sqrt(x**2+y**2+z**2)
	x_plane = sqrt(2d0)*x/(r_sphere+z)
	y_plane = sqrt(2d0)*y/(r_sphere+z)
 
end function



integer function matrix_rotation_to_top(x,y,z, rot)
! computes matrix rot for rotation of cube
! nearest to (x,y,z) vertex to (0,0,R)
! next nearest to (2sqrt2,0,3)

	implicit none
	
	real(8), intent(in) :: x,y,z
	real(8), dimension(1:3,1:3), intent(inout) :: rot

	rot(1:3,1) =(/  1, -2, -1 /) / sqrt(6d0)
	rot(1:3,2) =(/ -1, 0,  -1 /) / sqrt(2d0)
	rot(1:3,3) =(/  1, 1, -1 /) / sqrt(3d0)

	rot = transpose(rot)
	matrix_rotation_to_top = 1
end function matrix_rotation_to_top



integer function matrix_verge_rotation(x,y,z, rot)
	use matmul_module
	use simple_rotations
	implicit none

	interface
		 integer function index_verge_calc(x,y,z)
			 real(8) x,y,z
		 end function index_verge_calc
	end interface

	real(8) x,y,z
	real(8) rot(1:3,1:3)
	
	integer index_verge

	index_verge = index_verge_calc(x,y,z)
	matrix_verge_rotation = 0

	select case (index_verge)
		case(1)
			 rot=0d0; rot(1,1) = 1d0; rot(2,2)=1d0; rot(3,3)=1d0 ! identify matrix
		case(2)
			 call matmul1(transpose(Rotation_xy), Rotation_yz, rot)
		case(3)
			 rot = Rotation_yz
		case(4)
			 call matmul1(Rotation_xy, Rotation_yz, rot)
		case(5)
			 call matmul1(Rotation_xy, Rotation_xy, Rotation_yz, rot)
		case(6)
			 call matmul1(Rotation_yz, Rotation_yz, rot)
		case default
			 matrix_verge_rotation = 1
	endselect

	rot = transpose(rot)
end function matrix_verge_rotation
