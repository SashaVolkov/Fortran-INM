module morphism

interface
	 integer function cube2sphere(x, y, z, x_edge, y_edge, r_sphere, number_edge)
		 real(8) x, y, z
		 real(8) x_edge, y_edge, r_sphere
		 integer number_edge
	 end function cube2sphere

	 integer function reverse_stereo_projection(x_plane, y_plane, r_sphere, x, y, z)
		 real(8) x_plane, y_plane, r_sphere, x, y, z
	 end function reverse_stereo_projection

	 integer function stereo_projection(x_plane, y_plane, r_sphere, x, y, z)
		 real(8) x_plane, y_plane, r_sphere, x, y, z
	 end function stereo_projection
	 
	 integer function matrix_rotation_to_top(x,y,z, rot)
		 real(8) x,y,z
		 real(8), dimension(1:3,1:3) :: rot
	 end function matrix_rotation_to_top

	 integer function matrix_verge_rotation(x,y,z, rot)
		 real(8) x,y,z
		 real(8), dimension(1:3,1:3) :: rot
	 end function matrix_verge_rotation

end interface
end module morphism

module simple_rotations
	real(8), dimension(1:3,1:3) :: Rotation_xy = reshape(source = &
					 (/0, 1, 0, &
						-1, 0, 0, &
						 0, 0, 1 /), shape=(/3,3/))
	
	real(8), dimension(1:3,1:3) :: Rotation_yz = reshape(source = &
			 (/ (/ 1, 0, 0 /), &
				 (/ 0, 0, 1 /), &
				 (/ 0,-1, 0 /)  /), shape = (/3,3/))

	real(8), dimension(1:3,1:3) ::  Rotation_mir = reshape(source = &
			 (/ (/ 1, 0, 0 /), &
			 (/ 0,-1, 0 /), &
			 (/ 0, 0, 1 /)  /), shape = (/3,3/))
end module simple_rotations

module matrix_rotation
	interface
		 integer function init_compute_matrices(rots)
			 real(8) rots(3,3,48)
		 end function init_compute_matrices
	 integer function index_rotation_matrix(x,y,z)
		 real(8) x,y,z
	 end function index_rotation_matrix
	end interface
end module matrix_rotation

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


integer function matrix_vertex_rotation(x,y,z, rot)
	use matmul_module
	use simple_rotations
	implicit none

	interface
		 integer function index_vertex_calc(x,y,z)
			 real(8) x,y,z
		 end function index_vertex_calc
	end interface

	real(8) x,y,z
	real(8) rot(1:3,1:3)

	integer index_vertex

	index_vertex = index_vertex_calc(x,y,z)
	matrix_vertex_rotation = 0

	select case (index_vertex)
		 case(1)
				rot=0d0; rot(1,1)=1d0; rot(2,2)=1d0; rot(3,3)=1d0!identify matrix
		 case(2)
				call matmul1(Rotation_xy, Rotation_mir, rot)
		 case(3)
				rot = Rotation_xy
		 case(4)
				call matmul1(Rotation_xy, Rotation_xy, Rotation_mir, rot)
		 case(5)
				call matmul1(Rotation_xy, Rotation_xy, rot)
		 case(6)
				call matmul1(transpose(Rotation_xy),Rotation_mir, rot)
		 case(7)
				rot = transpose(Rotation_xy)
		 case(8)
				rot = Rotation_mir
		 case default
				matrix_vertex_rotation = 1
	end select
	rot=transpose(rot)
	matrix_vertex_rotation = 1

end function matrix_vertex_rotation


integer function index_vertex_calc(x,y,z)
	real(8) x,y,z

	index_vertex_calc = 0

	If ( x .ge.  y  .and.    y  .ge.  0d0) then
		 index_vertex_calc = 1
	else if( y .ge.  x  .and.    x  .ge.  0d0) then
		 index_vertex_calc = 2
	else if( y .ge. -x  .and.   -x  .ge.  0d0) then
		 index_vertex_calc = 3
	else if(-x .ge.  y  .and.    y  .ge.  0d0) then
		 index_vertex_calc = 4
	else if(-x .ge. -y  .and.   -y  .ge.  0d0) then
		 index_vertex_calc = 5
	else if(-y .ge.  x  .and.   -x  .ge.  0d0) then
		 index_vertex_calc = 6
	else if(-y .ge.  x  .and.    x  .ge.  0d0) then
		 index_vertex_calc = 7
	else if( x .ge. -y  .and.   -y  .ge.  0d0) then
		 index_vertex_calc = 8
	end if

end function index_vertex_calc

integer function index_verge_calc(x,y,z)

	real(8) x,y,z

	If ( (-z) .ge. max(abs(x),abs(y)) ) then
		 index_verge_calc = 1
	else if(   x  .ge. max(abs(y),abs(z)) ) then
		 index_verge_calc = 2
	else if(   y  .ge. max(abs(x),abs(z)) ) then
		 index_verge_calc = 3
	else if(  -x  .ge. max(abs(y),abs(z)) ) then
		 index_verge_calc = 4
	else if(  -y  .ge. max(abs(x),abs(z)) ) then
		 index_verge_calc = 5
	else if(   z  .ge. max(abs(x),abs(y)) ) then
		 index_verge_calc = 6
	end if

end function index_verge_calc

integer function index_rotation_matrix(x,y,z)
	interface
		 integer function index_verge_calc(x,y,z)
			 real(8) x,y,z
		 end function index_verge_calc

		 integer function index_vertex_calc(x,y,z)
			 real(8) x,y,z
		 end function index_vertex_calc

		 integer function matrix_verge_rotation(x,y,z,rot)
			 real(8) x,y,z, rot(1:3,1:3)
		 end function matrix_verge_rotation
	endinterface

	real(8) x,y,z
	real(8) rot(1:3,1:3), r(1:3)
	integer status

	index_rotation_matrix = (index_verge_calc(x,y,z)-1) * 8
	status = matrix_verge_rotation(x,y,z,rot)
	r = matmul(rot, (/x,y,z/))
	index_rotation_matrix = index_rotation_matrix + index_vertex_calc(r(1),r(2),r(3))
	

end function index_rotation_matrix


integer function init_compute_matrices(rots)
! computes matrices of rotation
	implicit none

	interface
		integer function cube2sphere(x,y,z,xe,ye,r,i)
			real(8) x,y,z,xe,ye,r
			integer i
		end function cube2sphere

		integer function matrix_verge_rotation(x,y,z,rot)
			real(8) x,y,z
			real(8) rot(3,3)
		end function matrix_verge_rotation

		integer function matrix_vertex_rotation(x,y,z,rot)
			real(8) x,y,z
			real(8) rot(3,3)
		end function matrix_vertex_rotation

		integer function matrix_rotation_to_top(x,y,z,rot)
			real(8) x,y,z
			real(8) rot(3,3)
		end function matrix_rotation_to_top

		integer function index_rotation_matrix(x,y,z)
			real(8) x,y,z
		end function index_rotation_matrix
	end interface


	real(8), dimension(1:3,1:3, 48) :: rots

	real(8), parameter :: pi2 = 6283185307179586d-15
	real(8) rot(3,3), r(3)
	real(8) x,y,z, x_edge, y_edge
	integer i, j, index, status
	
	rots=0
	do i = 1,6
		do j = 1,12
			if (mod(j,3).ne.0) then
				x_edge = cos(pi2*dble(j)/12d0)
				y_edge = sin(pi2*dble(j)/12d0)
				status = cube2sphere(x, y, z, x_edge,y_edge, 1d0, i)

				index = index_rotation_matrix(x,y,z)

				status = matrix_verge_rotation(x,y,z, rot)
				rots(1:3,1:3,index) = rot

				r = matmul(rot,(/x,y,z/))
				x=r(1);y=r(2);z=r(3)

				if ( matrix_vertex_rotation(x,y,z,rot) ) then
					rots(1:3,1:3,index) = matmul(rot, rots(1:3,1:3,index))
				end if 

				r = matmul(rot,(/x,y,z/))
				x=r(1); y=r(2); z=r(3)
				 
				if ( matrix_rotation_to_top(x,y,z,rot) ) then
					rots(1:3,1:3,index) = matmul(rot, rots(1:3,1:3,index))
				end if 

			end if
		end do
	end do
	init_compute_matrices = 1
end function init_compute_matrices
