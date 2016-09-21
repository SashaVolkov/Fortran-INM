module grid_interp

	use scan_print, Only: printer

implicit none

	Private
	Public :: interp

	Type interp

		Real(8), Allocatable :: latlon_c_off(:, :, :, :)
		Real(8), Allocatable :: latlon_c_to(:, :, :)
		Real(8), Allocatable :: weight(:, :, :)
		Real(8), Allocatable :: surface_off(:, :, :)
		Real(8), Allocatable :: surface_to(:, :)
		integer(4), Allocatable :: indexes_xyface(:, :, :, :)

		integer(4) dim, lon_max, lat_max

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Private :: alloc => alloc
		Procedure, Public :: deinit => deinit
	End Type


CONTAINS


	subroutine init(this, dim)

		Class(interp) :: this
		integer(4), intent(in) :: dim

		this.dim = dim;  this.lon_max = 360;  this.lat_max = 180

		call this.alloc()

	end subroutine



	subroutine alloc(this)
		Class(interp) :: this
		integer(4) dim, f, l, lon, lat

		dim = this.dim;  f = 1; l = 2*dim
		lon = this.lon_max; lat = this.lat_max

		Allocate(this.latlon_c_off(1:2, f:l, f:l, 1:6))
		Allocate(this.latlon_c_to(1:2, 1:lat, 1:lon))
		Allocate(this.surface_off(f:l, f:l, 1:6))
		Allocate(this.surface_to(1:lat, 1:lon))
		Allocate(this.weight(1:4, 1:lat, 1:lon))
		Allocate(this.indexes_xyface(1:3, 1:4, 1:lat, 1:lon))

	end subroutine



	subroutine deinit(this)
		Class(interp) :: this
		if (Allocated(this.latlon_c_off)) Deallocate(this.latlon_c_off)
		if (Allocated(this.latlon_c_to)) Deallocate(this.latlon_c_to)
		if (Allocated(this.surface_off)) Deallocate(this.surface_off)
		if (Allocated(this.surface_to)) Deallocate(this.surface_to)
		if (Allocated(this.weight)) Deallocate(this.weight)
		if (Allocated(this.indexes_xyface)) Deallocate(this.indexes_xyface)
	end subroutine



end module