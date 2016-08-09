module grid_generator_solver

use projections, Only: projection
use matrix_rotation, Only: matrix
use spherical
! use omp_lib

!			|6|
!		|5|2|3|4|    Baiburin diploma pages 5-8
!			|1|

implicit none

Private
Public :: generator

Type generator
	CONTAINS
	Procedure, Public :: conformal_cubed_sphere => conformal_cubed_sphere
	Procedure, Public :: equiangular_cubed_sphere => equiangular_cubed_sphere
	Procedure, Private :: face_to_octant => face_to_octant
		Procedure, Private :: rescale => rescale

End Type

CONTAINS

	subroutine conformal_cubed_sphere(this, dim, ns_xy, nf_xy, r_sphere, rescale, latlon_c, latlon)
		Class(generator) :: this
		integer, intent(in) :: dim, rescale, ns_xy(2), nf_xy(2)
		real(8), intent(in) :: r_sphere
		real(8), intent(out) :: latlon_c(1:2, 1:2*dim, 1:2*dim, 1:6)
		real(8), intent(out) :: latlon(1:2, 1:2*dim+1, 1:2*dim+1, 1:6)

		character*14 filename
		character istring
		integer(4) face_index, j, k, status, index, x_min, x_max, y_min, y_max, channel, l, m

		real(8) matr_of_rots(3,3,48), rot(1:3,1:3), r_vector(1:3),  t(2),  pi
		real(8) x_face,y_face, x_octant, y_octant, x_tan, y_tan, radius, longitude, latitude
		complex*16 w

		Type(projection) :: projection
		Type(matrix) :: matr

		pi = 314159265358979323846d-20


		call matr.compute_matr_of_rot(matr_of_rots, r_sphere) ! matr_of_rots - 48 2dim matrices (3*3). You can find them in grid/matrices_of_rotations.dat
		x_min = -2*dim; x_max = 2*dim; y_min = -2*dim; y_max = 2*dim

		do face_index= 1, 6
			do j= x_min, x_max
				do k= y_min, y_max

						x_face = j/dble(2*dim); y_face = k/dble(2*dim)		!normalized x coordinate

						if (rescale /= 0) then
							call this.rescale(x_face, y_face, x_face, y_face, rescale)
						end if

						call projection.stereographic_cube_to_sphere( r_vector, x_face, y_face, r_sphere, face_index, status) ! out :: r_vector, status !! Rancic p.978 (ii)
						call matr.index_rotation( r_vector, index) ! out :: index = from 1 to 48  !! 48 - full group of symetries of the cube

						call this.face_to_octant(x_face, y_face, x_octant, y_octant) ! first octant x > y > 0

						call projection.conformal_z_w( dcmplx( x_octant, y_octant), w)			! z = x + i*y !! Baiburin p.17 (3-5) !! Rancic p.978 (iii-iv)
						call projection.inverse( dreal(w), dimag(w), r_sphere, r_vector, status)	! out :: r_vector, status  !! Baiburin p.17 (6) !! Rancic p.978 (v)

						r_vector = matmul(transpose(matr_of_rots(1:3,1:3,index)),r_vector)		!! Baiburin p.17 (7) !! Rancic p.978 (v)
						! grid_points_xyz(:, j, k, face_index) = r_vector

						call cart2sphere(r_vector(1), r_vector(2), r_vector(3), radius, latitude, longitude)

						if(abs(mod(j,2)) == 1 .and. abs(mod(k,2)) == 1) then
							latlon_c(1, dim + (j+1)/2, dim + (k+1)/2, face_index) = latitude
							latlon_c(2, dim + (j+1)/2, dim + (k+1)/2, face_index) = longitude
						else if(abs(mod(j,2)) == 0 .and. abs(mod(k,2)) == 0) then
							latlon(1, dim + j/2 + 1, dim + k/2 + 1, face_index) = latitude
							latlon(2, dim + j/2 + 1, dim + k/2 + 1, face_index) = longitude
						end if

					end do
			end do
		end do


	end subroutine conformal_cubed_sphere



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



	subroutine rescale(this, x, y, x_tan, y_tan, rescale_type) ! x, y must be from 1 to -1
		Class(generator) :: this
		real(8), intent(in) :: x, y
		integer(4), intent(in) :: rescale_type
		real(8), intent(out) :: x_tan, y_tan

		if ( abs(x) > 1d0 .or. abs(y) > 1d0 ) then
			print *, x, y
		else if(rescale_type == 1) then
			x_tan = (1/dtan(2d0/3d0)) * dtan(x*2d0/3d0) ! DATAN like atan, but real(8)
			y_tan = (1/dtan(2d0/3d0)) * dtan(y*2d0/3d0)
		else if(rescale_type == 2) then
			x_tan = dsign(abs(x)**(65d0/51d0), x)
			y_tan = dsign(abs(y)**(65d0/51d0), y)
		end if



	end subroutine rescale




	subroutine equiangular_cubed_sphere(this, dim, step, equiang_c, latlon_c, latlon)

		Class(generator) :: this
		Type(projection) :: projection
		integer, intent(in) :: dim, step
		real(8), intent(out) :: latlon_c(1:2, 1:2*dim, 1:2*dim, 1:6)
		real(8), intent(out) :: equiang_c(1:2, 1:2*dim, 1:2*dim, 1:6)
		real(8), intent(out) :: latlon(1:2, 1:2*dim+1, 1:2*dim+1, 1:6)
		integer(4) face, i, j, k, min, max
		real(8) x,y,z,a, pi, r_vector(3), alpha,beta, latitude, longitude

		pi = 314159265358979323846d-20;  min = -2*dim;  max = - min


		do face = 1, 6
			do j= min, max
				do i= min, max
					alpha = pi*i/(8d0*dim);  beta= pi*j/(8d0*dim); k = sign(1, face - 3)

					latitude = datan(dtan(beta)*dcos(alpha));  longitude = alpha + (pi/2d0)*(face - 2)
					if ( face == 1 .or. face == 6 ) then
						latitude = -k*datan(dtan(alpha)/dtan(beta));  longitude = -datan(dcos(latitude)/dtan(beta))
					end if
						if(abs(mod(i,2)) == 1 .and. abs(mod(j,2)) == 1) then
							equiang_c(1, dim + (i+1)/2, dim + (j+1)/2, face) = alpha
							equiang_c(2, dim + (i+1)/2, dim + (j+1)/2, face) = beta
							latlon_c(1, dim + (i+1)/2, dim + (j+1)/2, face) = latitude
							latlon_c(2, dim + (i+1)/2, dim + (j+1)/2, face) = longitude
						else if(abs(mod(i,2)) == 0 .and. abs(mod(j,2)) == 0) then
							latlon(1, dim + i/2 + 1, dim + j/2 + 1, face) = latitude
							latlon(2, dim + i/2 + 1, dim + j/2 + 1, face) = longitude
						end if


						if ( i==0 .and. j==0 ) print *, latlon(:, dim + 1, dim+1, face), face

				end do
			end do
		end do


	end subroutine equiangular_cubed_sphere



end module