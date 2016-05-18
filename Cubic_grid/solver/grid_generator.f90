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
		Procedure, Private :: rescale => rescale

End Type

CONTAINS

	subroutine conformal_cubed_sphere_grid_generation(this, dim, r_sphere, rescale, grid_points_latlon)
		Class(generator) :: this
		integer, intent(in) :: dim, rescale
		real(8), intent(in) :: r_sphere
		real(8), intent(out) :: grid_points_latlon(1:2, -dim:dim, -dim:dim, 1:6) ! face_id, j, k, r_vector
		! real(8), intent(out) :: grid_points_xyz(1:3, -x_dimension:x_dimension, -y_dimension:y_dimension, 1:6) ! face_id, j, k, r_vector

		character*14 filename
		character istring
		integer(4) face_index, j, k, status, index, x_min, x_max, y_min, y_max, channel, l, m

		real(8) matr_of_rots(3,3,48), rot(1:3,1:3), r_vector(1:3),  t(2),  pi
		real(8) x_face,y_face, x_octant, y_octant, x_tan, y_tan, radius, longitude, latitude
		complex*16 w

		Type(projection) :: projection
		Type(matrix) :: matr

		pi = 314159265358979323846d-20

		call cpu_time(t(1)) ! Time start

		call matr.compute_matr_of_rot(matr_of_rots, r_sphere) ! matr_of_rots - 48 2dim matrices (3*3). You can find them in grid/matrices_of_rotations.dat
		x_min = -dim; x_max = dim; y_min = -dim; y_max = dim


		do face_index= 1, 6
			do j= x_min, x_max
! 				if ( abs(j) /=1 ) then
				do k= y_min, y_max
! 					if ( abs(k)/=1 ) then

						x_face = j/dble(dim); y_face = k/dble(dim)		!normalized x coordinate

						if (rescale == 1) then
							call this.rescale(x_face, y_face, x_face, y_face)
						end if

						call projection.stereographic_cube_to_sphere( r_vector, x_face, y_face, r_sphere, face_index, status) ! out :: r_vector, status !! Rancic p.978 (ii)
						call matr.index_rotation( r_vector, index) ! out :: index = from 1 to 48  !! 48 - full group of symetries of the cube

						call this.face_to_octant(x_face, y_face, x_octant, y_octant) ! first octant x > y > 0

						call projection.conformal_z_w( dcmplx( x_octant, y_octant), w)			! z = x + i*y !! Baiburin p.17 (3-5) !! Rancic p.978 (iii-iv)
						call projection.inverse( dreal(w), dimag(w), r_sphere, r_vector, status)	! out :: r_vector, status  !! Baiburin p.17 (6) !! Rancic p.978 (v)

						r_vector = matmul(transpose(matr_of_rots(1:3,1:3,index)),r_vector)		!! Baiburin p.17 (7) !! Rancic p.978 (v)
						! grid_points_xyz(:, j, k, face_index) = r_vector

						call cart2sphere(r_vector(1), r_vector(2), r_vector(3), radius, longitude, latitude)


						grid_points_latlon(1, j, k, face_index) = latitude
						grid_points_latlon(2, j, k, face_index) = longitude

! 						end if
					end do
! 				end if
			end do
		end do

		call cpu_time(t(2)) ! Time stop
! 		print '("")'
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



	subroutine rescale(this, x, y, x_tan, y_tan) ! x, y must be from 1 to -1
		Class(generator) :: this
		real(8), intent(in) :: x, y
		real(8), intent(out) :: x_tan, y_tan

		if ( abs(x) > 1d0 .or. abs(y) > 1d0 ) print *, x, y
		x_tan = (1/dtan(2d0/3d0)) * dtan(x*2d0/3d0) ! DATAN like atan, but real(8)
		y_tan = (1/dtan(2d0/3d0)) * dtan(y*2d0/3d0)
! 		x_tan = dsign(abs(x)**(4d0/3d0), x)
! 		y_tan = dsign(abs(y)**(4d0/3d0), y)


	end subroutine rescale





end module