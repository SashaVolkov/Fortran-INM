module interp_prec_to_cube


implicit none

	Private
	Public :: prec_to_cube

	Type prec_to_cube

		Real(8), Allocatable :: latlon_cubic(:, :, :, :)
		Real(8), Allocatable :: weight(:, :, :, :)
		Real(4), Allocatable :: surface_to(:, :, :)
		integer(4), Allocatable :: indexes_ll(:, :, :, :, :)

		integer(4) dim, lon_max, lat_max, step

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Private :: alloc => alloc
		Procedure, Public :: deinit => deinit
		Procedure, Public :: weight_find => weight_find
		Procedure, Public :: cell_search => cell_search
		Procedure, Public :: interpolate => interpolate
	End Type


CONTAINS


	subroutine init(this, dim, lon_max, lat_max)

		Class(prec_to_cube) :: this
		integer(4), intent(in) :: dim, lon_max, lat_max

		this.dim = dim;  this.lon_max = lon_max;  this.lat_max = lat_max
		this.step = 2

		call this.alloc()

	end subroutine



	subroutine alloc(this)
		Class(prec_to_cube) :: this
		integer(4) dim, f, l, lon, lat

		dim = this.dim;  f = 1 - this.step; l = 2*dim + this.step
		lon = this.lon_max; lat = this.lat_max

		Allocate(this.latlon_cubic(2, f:l, f:l, 6))
		Allocate(this.surface_to(1:2*dim, 1:2*dim, 6))
		Allocate(this.weight(4, 1:2*dim, 1:2*dim, 6))
		Allocate(this.indexes_ll(2, 4, 1:2*dim, 1:2*dim, 6))

	end subroutine



	subroutine deinit(this)
		Class(prec_to_cube) :: this
		if (Allocated(this.latlon_cubic)) Deallocate(this.latlon_cubic)
		if (Allocated(this.surface_to)) Deallocate(this.surface_to)
		if (Allocated(this.weight)) Deallocate(this.weight)
		if (Allocated(this.indexes_ll)) Deallocate(this.indexes_ll)
	end subroutine



	subroutine weight_find(this)
		Class(prec_to_cube) :: this
		integer(4) lon, lat, x, y, face, i, j, dim
		Real(8) :: S(4), Big_S, latlon(2), latlon1(2), latlon2(2), d(4), ang(4), d_xy(2,4), lagr(2,4), factor
		real(8), parameter :: pi = 314159265358979323846d-20

		call this.cell_search()

		dim = this.dim;  factor = 90d0/dble(this.lat_max)/180d0*pi

		do face = 1, 6
			do x = 1, 2*dim
				do y = 1, 2*dim

					d(1) = abs(this.latlon_cubic(1, x, y, face)/factor - this.indexes_ll(1,1, x, y, face))
					d(2) = abs(this.latlon_cubic(2, x, y, face)/factor - this.indexes_ll(2,2, x, y, face))
					d(3) = abs(this.latlon_cubic(1, x, y, face)/factor - this.indexes_ll(1,3, x, y, face))
					d(4) = abs(this.latlon_cubic(2, x, y, face)/factor - this.indexes_ll(2,4, x, y, face))

					Big_S = (d(1) + d(3))*(d(2) + d(4))

					this.weight(1, x, y, face) = d(2)*d(3)/Big_S
					this.weight(2, x, y, face) = d(4)*d(3)/Big_S
					this.weight(3, x, y, face) = d(1)*d(4)/Big_S
					this.weight(4, x, y, face) = d(2)*d(1)/Big_S

				end do
			end do
		end do


		! print *, d(:), Big_S
	end subroutine




	subroutine interpolate(this, surface_off)
		Class(prec_to_cube) :: this
		integer(4) lon(4), lat(4), x, y, face
		Real(4), Intent(in) :: surface_off(-this.lon_max:this.lon_max, -this.lat_max:this.lat_max)
		Real(8) :: weight(4)

		do face = 1, 6
			do x = 1, 2*this.dim
				do y = 1, 2*this.dim

				lat(:) = this.indexes_ll(1,:, x, y, face)
				lon(:) = this.indexes_ll(2,:, x, y, face)
				weight(:) = this.weight(:, x, y, face)

this.surface_to(x, y, face) = surface_off(lon(1), lat(1))*weight(1) + surface_off(lon(2), lat(2))*weight(2) + surface_off(lon(3), lat(3))*weight(3) + surface_off(lon(4), lat(4))*weight(4)

				end do
			end do
		end do

	end subroutine




	subroutine cell_search(this)
		Class(prec_to_cube) :: this
		integer(4) dim, x, y, face
		Real(8) :: factor, f, lat_cell(2), lon_cell(2)
		real(8), parameter :: pi = 314159265358979323846d-20
		dim = this.dim;  f = 90d0/dble(this.lat_max);  factor = f/180d0*pi

		do face = 1, 6
			do x = 1, 2*dim
				do y = 1, 2*dim
					lat_cell(1) = ceiling(this.latlon_cubic(1, x, y, face)/factor)
					lat_cell(2) = floor(this.latlon_cubic(1, x, y, face)/factor)
					lon_cell(1) = ceiling(this.latlon_cubic(2, x, y, face)/factor)
					lon_cell(2) = floor(this.latlon_cubic(2, x, y, face)/factor)
					if(lon_cell(2) == lon_cell(1)) lon_cell(2) = lon_cell(1) - 1
					if(lon_cell(2) == -180/f) lon_cell(2) = - lon_cell(2)


					this.indexes_ll(1,1, x, y, face) = lat_cell(1);  this.indexes_ll(2,1, x, y, face) = lon_cell(2)
					this.indexes_ll(1,2, x, y, face) = lat_cell(1);  this.indexes_ll(2,2, x, y, face) = lon_cell(1)
					this.indexes_ll(1,3, x, y, face) = lat_cell(2);  this.indexes_ll(2,3, x, y, face) = lon_cell(1)
					this.indexes_ll(1,4, x, y, face) = lat_cell(2);  this.indexes_ll(2,4, x, y, face) = lon_cell(2)

				end do
			end do
		end do

	end subroutine



end module