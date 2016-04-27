module grid_generator_solver

use projections, Only: projection
use matrix_rotation, Only: matrix
use spherical

!			|6|
!		|5|2|3|4|    Baiburin diploma pages 5-8
!			|1|

implicit none

Private
Public :: generator

Type generator
	CONTAINS
	Procedure, Public :: conformal_cubed_sphere => conformal_cubed_sphere_grid_generation
	Procedure, Private :: face_to_octant => face_to_octant
	Procedure, Private :: Adcroft_tan => Adcroft_tan
End Type

CONTAINS

	subroutine conformal_cubed_sphere_grid_generation(this, x_dimension,y_dimension, r_sphere, grid_points)
		Class(generator) :: this
		integer, intent(in) :: x_dimension, y_dimension
		real(8), intent(in) :: r_sphere
		real(8), intent(out) :: grid_points(1:2, -x_dimension:x_dimension, -y_dimension:y_dimension, 1:6) ! face_id, j, k, r_vector

		character*14 filename
		character istring
		integer face_index, j, k, status, index, x_min, x_max, y_min, y_max, channel

		real(8) matr_of_rots(3,3,48), rot(1:3,1:3), r_vector(1:3),  t(2)
		real(8) x_face,y_face, x_octant, y_octant, x_tan, y_tan, radius, longitude, latitude
		complex*16 w

		Type(projection) :: projection
		Type(matrix) :: matr

		call cpu_time(t(1)) ! Time start

		call matr.compute_matr_of_rot(matr_of_rots, r_sphere) ! matr_of_rots - 48 2dim matrices (3*3). You can find them in grid/matrices_of_rotations.dat
		x_min = -x_dimension; x_max = x_dimension; y_min = -y_dimension; y_max = y_dimension


		do face_index= 1, 6
			do j= x_min, x_max
				do k= y_min, y_max

					x_face = j/dble(x_dimension); y_face = k/dble(y_dimension)		!normalized x coordinate

					call projection.stereographic_cube_to_sphere( r_vector, x_face, y_face, r_sphere, face_index, status) ! out :: r_vector, status !! Rancic p.978 (ii)
					call matr.index_rotation( r_vector, index) ! out :: index = from 1 to 48  !! 48 - full group of symetries of the cube

					call this.face_to_octant(x_face, y_face, x_octant, y_octant) ! first octant x > y > 0

					call projection.conformal_z_w( dcmplx( x_octant, y_octant), w)			! z = x + i*y !! Baiburin p.17 (3-5) !! Rancic p.978 (iii-iv)
					call projection.inverse( dreal(w), dimag(w), r_sphere, r_vector, status)	! out :: r_vector, status  !! Baiburin p.17 (6) !! Rancic p.978 (v)

					r_vector = matmul(transpose(matr_of_rots(1:3,1:3,index)),r_vector)		!! Baiburin p.17 (7) !! Rancic p.978 (v)

					call cart2sphere(r_vector(1), r_vector(2), r_vector(3), radius, longitude, latitude)

						grid_points(1, j, k, face_index) = latitude
						grid_points(2, j, k, face_index) = longitude

				end do
			end do
		end do

		call cpu_time(t(2)) ! Time stop
		print '("")'
		print '(" Time of generation = ", f6.3, " sec")', t(2) - t(1)

	end subroutine conformal_cubed_sphere_grid_generation



	subroutine face_to_octant(this, x_face, y_face, x_octant, y_octant)
	!x_octant >= 0; y_octant >= 0 and x_octant > y_octant - first octant
		Class(generator) :: this
		real(8), intent(out) :: x_octant, y_octant
		real(8), intent(in) :: x_face, y_face

					if (abs(x_face)>abs(y_face)) then
						y_octant = abs(1-abs(x_face))/2d0
						x_octant = abs(1-abs(y_face))/2d0
					else
						x_octant = abs(1-abs(x_face))/2d0				! x = |(1 - |x|)|/2
						y_octant = abs(1-abs(y_face))/2d0				! y = |(1 - |y|)|/2
					end if

	end subroutine face_to_octant



	subroutine Adcroft_tan(this, x, y, x_tan, y_tan) ! x, y must be from 1 to -1
		Class(generator) :: this
		real(8), intent(in) :: x, y
		real(8), intent(out) :: x_tan, y_tan

		x_tan = (1/datan(2d0/3d0)) * datan(x*2d0/3d0) ! DATAN like atan, but real(8)
		y_tan = (1/datan(2d0/3d0)) * datan(y*2d0/3d0)

	end subroutine Adcroft_tan



end module