module matrix_rotation

use projections, Only: projection
use matmul_module
use simple_rotations

implicit none

	Private
	Public :: matrix

	Type matrix
		CONTAINS
		Procedure, Public :: compute_matr_of_rot => compute_matr_of_rot
		Procedure, Public :: index_rotation => index_rotation_matrix
	End Type


CONTAINS

	subroutine compute_matr_of_rot(this, matr_of_rots)
	! matr_of_rots - 48 2dim matrices (3*3)

		real(8), intent(out) :: matr_of_rots(1:3,1:3, 48)

		real(8), parameter :: pi2 = 6283185307179586d-15
		real(8) rot(3,3), r_vector(1:3), r_sphere
		real(8) x_face, y_face
		integer i, j, index, status
		Type(projection) :: projection
		Class(matrix) :: this

		r_sphere = 1d0
		
		matr_of_rots=0
		do i = 1,6
			do j = 1,12
				if (mod(j,3).ne.0) then						!MOD(x,y) = remainder x - INT(x/y)*y
					x_face = cos(pi2*dble(j)/12d0)	!DBLE(A) Converts A to double precision real type
					y_face = sin(pi2*dble(j)/12d0)
					call projection.stereographic_cube_to_sphere( r_vector, x_face,y_face, r_sphere, i, status)

					call this.index_rotation(r_vector,index)	! Calculating index. It can be from 1 to 48.

					call matrix_face_rotation(r_vector, rot)
					matr_of_rots(1:3,1:3,index) = rot

					r_vector = matmul(rot,r_vector)

					call matrix_octant_rotation(r_vector,rot)
					matr_of_rots(1:3,1:3,index) = matmul(rot, matr_of_rots(1:3,1:3,index))

					r_vector = matmul(rot,r_vector)

					call matrix_rotation_to_top(rot)
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

	end subroutine compute_matr_of_rot



	subroutine index_rotation_matrix(this, r_vector, index)
	! Calculating index. It can be from 1 to 48.

		Class(matrix) :: this
		real(8), intent(in) :: r_vector(1:3)
		real(8) rot(1:3,1:3), r(1:3)
		integer(4), intent(out) :: index

		index = (select_face_of_cube_index(r_vector)-1) * 8			! 0*8, 1*8, ... , 5*8 - one of them
		call matrix_face_rotation(r_vector,rot)
		r = matmul(rot, r_vector)
		index = index + select_octant_of_face_index(r)		! index + (1 or 2 or ... or 8) ! max index = 48, min index = 1

	end subroutine index_rotation_matrix



	subroutine matrix_rotation_to_top(rot)
	! computes matrix rot for rotation of cube
	! nearest to (x,y,z) vertex to (0,0,R)
	! next nearest to (2sqrt2,0,3)

		real(8), intent(out) :: rot(1:3,1:3)

		rot(1:3,1) =(/  1, -2, -1 /) / sqrt(6d0)
		rot(1:3,2) =(/ -1, 0,  -1 /) / sqrt(2d0)
		rot(1:3,3) =(/  1, 1, -1 /) / sqrt(3d0)

		rot = transpose(rot)

	end subroutine matrix_rotation_to_top



	subroutine matrix_face_rotation(r_vector, rot)

		real(8), intent(out) :: rot(1:3,1:3)
		real(8), intent(in) :: r_vector(1:3)
		integer face_index

		face_index = select_face_of_cube_index(r_vector)								! Define one of 6 faces of the cube

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
		end select

		rot = transpose(rot)

	end subroutine matrix_face_rotation



	subroutine matrix_octant_rotation(r_vector,rot)

		real(8), intent(out) :: rot(1:3,1:3)
		real(8), intent(in) :: r_vector(1:3)
		integer octant_index

		octant_index = select_octant_of_face_index(r_vector)

		select case (octant_index)
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
		end select

		rot=transpose(rot)

	end subroutine matrix_octant_rotation



	integer function select_face_of_cube_index(r_vector) !Look for "Net of cube"

		real(8), intent(in) :: r_vector(1:3)
		real(8) x,y,z
		x = r_vector(1); y = r_vector(2); z = r_vector(3)

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



	integer function select_octant_of_face_index(r_vector)

		real(8), intent(in) :: r_vector(1:3)
		real(8) x,y,z
		x = r_vector(1); y = r_vector(2); z = r_vector(3)

		select_octant_of_face_index = 0

		If ( x >=  y  .and.    y  >=  0d0) then
			 select_octant_of_face_index = 1
		else if( y >=  x  .and.    x  >=  0d0) then
			 select_octant_of_face_index = 2
		else if( y >= -x  .and.   -x  >=  0d0) then
			 select_octant_of_face_index = 3
		else if(-x >=  y  .and.    y  >=  0d0) then
			 select_octant_of_face_index = 4
		else if(-x >= -y  .and.   -y  >=  0d0) then
			 select_octant_of_face_index = 5
		else if(-y >=  x  .and.   -x  >=  0d0) then
			 select_octant_of_face_index = 6
		else if(-y >=  x  .and.    x  >=  0d0) then
			 select_octant_of_face_index = 7
		else if( x >= -y  .and.   -y  >=  0d0) then
			 select_octant_of_face_index = 8
		end if

	end function select_octant_of_face_index

end module