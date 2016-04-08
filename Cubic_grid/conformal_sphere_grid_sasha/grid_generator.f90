module grid_generator

use projections, Only: projection
use matrix_rotation, Only: matrix

implicit none

Private
Public :: grid

Type grid
	CONTAINS
	Procedure, Public :: conformal_cubed_sphere => conformal_cubed_sphere_grid_generation
	Procedure, Private :: face_to_triangle => face_to_triangle
End Type

CONTAINS


	subroutine conformal_cubed_sphere_grid_generation(this, x_dimension,y_dimension)
		Class(grid) :: this
		integer, intent(in) :: x_dimension, y_dimension

		character*14 filename
		character istring
		integer face_index, j, k, status, index, x_min, x_max, y_min, y_max

		real(8) matr_of_rots(3,3,48), rot(1:3,1:3)
		real(8) r(1:3)
		real(8) x,y,z
		real(8) x_face,y_face, x_triangle, y_triangle
		real(8) r_sphere
		complex*16 w


		Type(projection) :: projection
		Type(matrix) :: matr

		call matr.init_matrices(matr_of_rots) ! matr_of_rots - 48 2dim matrices (3*3). You can find them in grid/matrices_of_rotations.dat
		x_min = -x_dimension; x_max = x_dimension; y_min = -y_dimension; y_max = y_dimension
		r_sphere = 1d0


		do face_index= 1, 6

			write(istring(1:1), '(i1.1)') face_index
			filename = "grid/face" // istring // ".dat" 
			open (20, file = filename)

			do j= x_min, x_max
				do k= y_min, y_max

					x_face = j/dble(x_dimension)		!normalized x coordinate
					y_face = k/dble(y_dimension)

					call projection.stereographic_cube_to_sphere( x, y, z, x_face, y_face, r_sphere, face_index, status) ! out :: x, y, z, status !! Rancic p.978 (ii)
					call matr.index_rotation( x, y, z, index) ! out :: index = from 1 to 48  !! 48 - full group of symetries of the cube

					call this.face_to_triangle(x_face, y_face, x_triangle, y_triangle)

					call projection.conformal_z_w( dcmplx( x_triangle, y_triangle), w)			! z = x_triangle + i*y_triangle  !! Baiburin p.17 (3-5) !! Rancic p.978 (iii-iv)
					call projection.inverse( dreal(w), dimag(w), r_sphere, x, y, z, status)	! intent(out) :: x, y, z, status  !! Baiburin p.17 (6)

					r = matmul(transpose(matr_of_rots(1:3,1:3,index)),(/x,y,z/))		!! Baiburin p.17 (7)
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



	subroutine face_to_triangle(this, x_face, y_face, x_triangle, y_triangle)
	!x_triangle >= 0; y_triangle >= 0 and x_triangle > y_triangle
	! это гарантии последнего предложения на 13й странице, а именно x<y , координаты на грани должны попасть в 1/8 часть, если нет, то мы преобразуем их по симметрии относительно прямой x=y.
		Class(grid) :: this
		real(8), intent(out) :: x_triangle, y_triangle
		real(8), intent(in) :: x_face, y_face

					if (abs(x_face)>abs(y_face)) then
						y_triangle = abs(1-abs(x_face))/2d0
						x_triangle = abs(1-abs(y_face))/2d0
					else
						x_triangle = abs(1-abs(x_face))/2d0				! x = |(1 - |x|)|/2
						y_triangle = abs(1-abs(y_face))/2d0				! y = |(1 - |y|)|/2
					end if

	end subroutine face_to_triangle

end module