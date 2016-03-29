module grid_generator

use projections, Only: projection
use matrix_rotation, Only: matrix

implicit none

Private
Public :: grid

Type grid
	CONTAINS
	Procedure :: conformal_cubed_sphere => conformal_cubed_sphere_grid_generation
End Type

CONTAINS


	subroutine conformal_cubed_sphere_grid_generation(this, x_dimension,y_dimension)
		Class(grid) :: this
		integer, intent(in) :: x_dimension, y_dimension
		character*14 filename
		character istring

		real(8) matr_of_rots(3,3,48)
		real(8) rot(1:3,1:3)
		real(8) r(1:3)
		real(8) x,y,z
		real(8) x_face,y_face, x_edge, y_edge
		real(8) r_sphere
		complex*16 w

		integer face_index, j, k, status, index, x_min, x_max, y_min, y_max
		Type(projection) :: projection
		Type(matrix) :: matr

		call matr.init_matrices(matr_of_rots) ! matr_of_rots - 48 2dim matrices (3*3). You can find them in grid/matrices_of_rotations.dat

		x_min = -x_dimension; x_max = x_dimension; y_min = -y_dimension; y_max = y_dimension

		do face_index= 1, 6

			write(istring(1:1), '(i1.1)') face_index
			filename = "grid/face" // istring // ".dat" 
			open (20, file = filename)

			do j= x_min, x_max
				do k= y_min, y_max

					x_face = j/dble(x_dimension)		!normalized x coordinate
					y_face = k/dble(y_dimension)
					r_sphere = 1d0

					call projection.cube2sphere( x, y, z, x_face, y_face, r_sphere, face_index, status) ! out :: x, y, z, status
					call matr.index_rotation( x, y, z, index) ! out :: index = from 1 to 48

					if (abs(x_face)>abs(y_face)) then
						y_edge = abs(sign(1d0,x_face)*(1-abs(x_face)))/2d0				! SIGN(A,B) returns the value of A with the sign of B.
						x_edge = abs(sign(1d0,y_face)*(1-abs(y_face)))/2d0				! WTF?
					else
						x_edge = abs(sign(1d0,x_face)*(1-abs(x_face)))/2d0				! sgn(x)*(1 - |x|)/2
						y_edge = abs(sign(1d0,y_face)*(1-abs(y_face)))/2d0
					end if

					! z = x_edge + i*y_edge

					call projection.conformal_z_w(dcmplx(x_edge,y_edge), w) 
! DCMPLX(X [,Y]) returns a double complex number where X is converted to the real component.
! If Y is present it is converted to the imaginary component if not present then imaginary is set to 0.0.
! If X is complex then Y must not be present. 
					call projection.inverse(dreal(w),dimag(w),r_sphere,x,y,z,status)

					r = matmul(transpose(matr_of_rots(1:3,1:3,index)),(/x,y,z/))
					x = r(1); y=r(2); z=r(3)

					write(20,*) x,y,z

				end do
				write(20,*)
			end do

		end do
		close(20)

		open(20, file = "grid/parameters.dat")
		write(20,*) x_dimension, y_dimension
		close(20)
	end subroutine conformal_cubed_sphere_grid_generation


end module