module matrix_rotation

use mapping, Only: mapp
use matmul_module
use simple_rotations

implicit none

	Private
	Public :: matrix

	Type matrix
		CONTAINS
		Procedure, Public :: init_matrices => init_compute_matrices
		Procedure, Private :: matrix_rotation_to_top => matrix_rotation_to_top
		Procedure, Public :: index_rotation => index_rotation_matrix
	End Type


CONTAINS

	subroutine init_compute_matrices(this, matr_of_rots)
	! computes 48 matrices of rotation
	! matr_of_rots - 48 2dim matrices (3*3)

		real(8), intent(out), dimension(1:3,1:3, 48) :: matr_of_rots

		real(8), parameter :: pi2 = 6283185307179586d-15
		real(8) rot(3,3), r(3)
		real(8) x,y,z, x_face, y_face
		integer i, j, index, status
		Type(mapp) :: map
		Class(matrix) :: this
		
		matr_of_rots=0
		do i = 1,6
			do j = 1,12
				if (mod(j,3).ne.0) then						!MOD(x,y) = remainder x - INT(x/y)*y
					x_face = cos(pi2*dble(j)/12d0)	!DBLE(A) Converts A to double precision real type
					y_face = sin(pi2*dble(j)/12d0)
					call map.cube2sphere(x, y, z, x_face,y_face, 1d0, i, status)

					call this.index_rotation(x,y,z,index)	! Calculating index. It can be from 1 to 48.

					status = matrix_face_rotation(x,y,z, rot)
					matr_of_rots(1:3,1:3,index) = rot

					r = matmul(rot,(/x,y,z/))
					x=r(1);y=r(2);z=r(3)

					status = matrix_vertex_rotation(x,y,z,rot)
						matr_of_rots(1:3,1:3,index) = matmul(rot, matr_of_rots(1:3,1:3,index))

					r = matmul(rot,(/x,y,z/))
					x=r(1); y=r(2); z=r(3)

					call this.matrix_rotation_to_top(x,y,z,rot,status)
					matr_of_rots(1:3,1:3,index) = matmul(rot, matr_of_rots(1:3,1:3,index))

				end if
			end do
		end do


		open (20, file = "grid/matrices_of_rotations.dat")
		do i =1,48
			 write(20,*) matr_of_rots(1:3,1,i)
			 write(20,*) matr_of_rots(1:3,2,i)
			 write(20,*) matr_of_rots(1:3,3,i)
			 write(20,*)
		end do
		close(20)

	end subroutine init_compute_matrices



	subroutine matrix_rotation_to_top(this, x ,y ,z ,rot, status)
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



	subroutine index_rotation_matrix(this, x, y, z, index)
	! Calculating index. It can be from 1 to 48.

		Class(matrix) :: this
		real(8), intent(in) :: x,y,z
		real(8) rot(1:3,1:3), r(1:3)
		integer(4), intent(out) :: index
		integer(4) status

		index = (select_face_of_cube_index(x,y,z)-1) * 8			! 0*8, 1*8, ... , 5*8 - one of them
		status = matrix_face_rotation(x,y,z,rot)
		r = matmul(rot, (/x,y,z/))
		index = index + select_vertex_index(r(1),r(2),r(3))		! index + (1 or 2 or ... or 8) ! max index = 48, min index = 1

	end subroutine index_rotation_matrix



	integer function matrix_face_rotation(x, y, z, rot)

		real(8) x,y,z
		real(8) rot(1:3,1:3)
		integer face_index

		face_index = select_face_of_cube_index(x,y,z)								! Define one of 6 faces of the cube
		matrix_face_rotation = 0

		select case (face_index)
			case(1)
				 rot=0d0; rot(1,1) = 1d0; rot(2,2)=1d0; rot(3,3)=1d0 			! identify matrix
			case(2)
				 call matmul1(transpose(Rotation_xy), Rotation_yz, rot)		! All Rotation_... defined in simple_rotations module
			case(3)
				 rot = Rotation_yz
			case(4)
				 call matmul1(Rotation_xy, Rotation_yz, rot)
			case(5)
				 call matmul1(Rotation_xy, Rotation_xy, Rotation_yz, rot)
			case(6)
				 call matmul1(Rotation_yz, Rotation_yz, rot)
			case default
				 matrix_face_rotation = 1
		end select

		rot = transpose(rot)
	end function matrix_face_rotation



	integer function matrix_vertex_rotation(x, y, z,rot)

		real(8) x,y,z
		real(8) rot(1:3,1:3)

		integer index_vertex

		index_vertex = select_vertex_index(x,y,z)
		matrix_vertex_rotation = 0

		select case (index_vertex)
			 case(1)
					rot=0d0; rot(1,1)=1d0; rot(2,2)=1d0; rot(3,3)=1d0 !identify matrix
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



	integer function select_face_of_cube_index(x, y, z) !Look for "Net of cube"

		real(8) x,y,z

!			|6|
!		|5|2|3|4|    Baiburin diploma pages 5-8
!			|1|

		If ( (-z) >= max(abs(x),abs(y)) ) then
			 select_face_of_cube_index = 1
		else if(   x  >= max(abs(y),abs(z)) ) then
			 select_face_of_cube_index = 2
		else if(   y  >= max(abs(x),abs(z)) ) then
			 select_face_of_cube_index = 3
		else if(  -x  >= max(abs(y),abs(z)) ) then
			 select_face_of_cube_index = 4
		else if(  -y  >= max(abs(x),abs(z)) ) then
			 select_face_of_cube_index = 5
		else if(   z  >= max(abs(x),abs(y)) ) then
			 select_face_of_cube_index = 6
		end if

	end function select_face_of_cube_index



	integer function select_vertex_index(x, y, z)
		real(8) x,y,z

		select_vertex_index = 0

		If ( x >=  y  .and.    y  >=  0d0) then
			 select_vertex_index = 1
		else if( y >=  x  .and.    x  >=  0d0) then
			 select_vertex_index = 2
		else if( y >= -x  .and.   -x  >=  0d0) then
			 select_vertex_index = 3
		else if(-x >=  y  .and.    y  >=  0d0) then
			 select_vertex_index = 4
		else if(-x >= -y  .and.   -y  >=  0d0) then
			 select_vertex_index = 5
		else if(-y >=  x  .and.   -x  >=  0d0) then
			 select_vertex_index = 6
		else if(-y >=  x  .and.    x  >=  0d0) then
			 select_vertex_index = 7
		else if( x >= -y  .and.   -y  >=  0d0) then
			 select_vertex_index = 8
		end if

	end function select_vertex_index

end module