module matrix_rotation

! use subfunc, Only: func
use mapping, Only: conf
use matmul_module
use simple_rotations

implicit none

	Public :: matrix

	Type matrix
		CONTAINS
		Procedure :: init_compute_matrices => init_compute_matrices
		Procedure :: matrix_rotation_to_top => matrix_rotation_to_top
	End Type

CONTAINS



	subroutine init_compute_matrices(this, rots)
	! computes matrices of rotation

		real(8), dimension(1:3,1:3, 48) :: rots

		real(8), parameter :: pi2 = 6283185307179586d-15
		real(8) rot(3,3), r(3)
		real(8) x,y,z, x_edge, y_edge
		integer i, j, index, status
		Type(conf) :: conformal
		Class(matrix) :: this
		
		rots=0
		do i = 1,6
			do j = 1,12
				if (mod(j,3).ne.0) then
					x_edge = cos(pi2*dble(j)/12d0)
					y_edge = sin(pi2*dble(j)/12d0)
					call conformal.cube2sphere(x, y, z, x_edge,y_edge, 1d0, i, status)

					index = index_rotation_matrix(x,y,z)

					status = matrix_verge_rotation(x,y,z, rot)
					rots(1:3,1:3,index) = rot

					r = matmul(rot,(/x,y,z/))
					x=r(1);y=r(2);z=r(3)

					status = matrix_vertex_rotation(x,y,z,rot)
						rots(1:3,1:3,index) = matmul(rot, rots(1:3,1:3,index))

					r = matmul(rot,(/x,y,z/))
					x=r(1); y=r(2); z=r(3)
					

					call this.matrix_rotation_to_top(x,y,z,rot,status)
					rots(1:3,1:3,index) = matmul(rot, rots(1:3,1:3,index))

				end if
			end do
		end do
	end subroutine init_compute_matrices

	subroutine matrix_rotation_to_top(this,x,y,z, rot, status)
	! computes matrix rot for rotation of cube
	! nearest to (x,y,z) vertex to (0,0,R)
	! next nearest to (2sqrt2,0,3)
		
		real(8), intent(in) :: x,y,z
		real(8), dimension(1:3,1:3), intent(out) :: rot
		integer(4), intent(out) :: status
		Class(matrix) :: this

		rot(1:3,1) =(/  1, -2, -1 /) / sqrt(6d0)
		rot(1:3,2) =(/ -1, 0,  -1 /) / sqrt(2d0)
		rot(1:3,3) =(/  1, 1, -1 /) / sqrt(3d0)

		rot = transpose(rot)
		status = 1
	end subroutine matrix_rotation_to_top



	integer function index_rotation_matrix(x,y,z)

		real(8) x,y,z
		real(8) rot(1:3,1:3), r(1:3)
		integer status

		index_rotation_matrix = (index_verge_calc(x,y,z)-1) * 8
		status = matrix_verge_rotation(x,y,z,rot)
		r = matmul(rot, (/x,y,z/))
		index_rotation_matrix = index_rotation_matrix + index_vertex_calc(r(1),r(2),r(3))

	end function index_rotation_matrix



	integer function matrix_verge_rotation(x,y,z, rot)

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
		end select

		rot = transpose(rot)
	end function matrix_verge_rotation



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




	integer function matrix_vertex_rotation(x,y,z, rot)

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


end module matrix_rotation