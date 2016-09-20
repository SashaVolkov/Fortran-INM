module grid_interp


implicit none

	Private
	Public :: interp

	Type interp

		Real(8), Allocatable :: latlon_c_off(:, :, :, :)
		Real(8), Allocatable :: latlon_c_to(:, :, :)
		Real(8), Allocatable :: weight(:, :, :)
		integer(4), Allocatable :: indexes_xyface(:, :, :)

		Real(8) ::  omega_cor, r_sphere, g, dt, dx_min, dy_min, dx_max, dy_max, pi, delta_on_cube
		integer(4) dim, step, rescale, ns_xy(2), nf_xy(2), grid_type
		integer(4) first_x, first_y, last_x, last_y, snd_xy(6, 4, 2), rcv_xy(6, 4, 2)

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Private :: alloc => alloc
		Procedure, Public :: deinit => deinit
		Procedure, Private :: const_def => const_def
	End Type


CONTAINS


	subroutine init(this, grid)

		Class(interp) :: this
		Class(g_var) :: grid

		this.dim = grid.dim;  this.step = grid.step
		this.ns_xy(:) = grid.ns_xy(:);  this.nf_xy(:) = grid.nf_xy(:)

		call this.alloc()

	end subroutine



	subroutine alloc(this)
		Class(interp) :: this
		integer(4) dim, step, f, l

		dim = this.dim;  step = this.step;  f = 1-step; l = 2*dim + step

		Allocate(this.latlon_c_off(1:2, f:l, f:l, 1:6))
		Allocate(this.latlon_c_to(1:2, 0:179, 0:359))
		Allocate(this.weight(1:4, 0:179, 0:359))
		Allocate(this.indexes_xyface(1:3, 0:179, 0:359))

	end subroutine



	subroutine deinit(this)
		Class(g_var) :: this
		if (Allocated(this.f_cor)) Deallocate(this.f_cor)
		if (Allocated(this.latlon_c)) Deallocate(this.latlon_c)
		if (Allocated(this.cube_coord_c)) Deallocate(this.cube_coord_c)
		if (Allocated(this.latlon)) Deallocate(this.latlon)

		if (Allocated(this.square)) Deallocate(this.square)
		if (Allocated(this.triangle_area)) Deallocate(this.triangle_area)
		if (Allocated(this.triangle_angles)) Deallocate(this.triangle_angles)
		if (Allocated(this.square_angles)) Deallocate(this.square_angles)
	end subroutine



	subroutine const_def(this, g, metr)
		Class(g_var) :: this
		Class(geometry) :: g
		Class(metric) :: metr
		real(8) dist, omega_cor, h(-1:2)
		integer(4) face, x, y, dim
		integer(4), parameter :: A =1, B=2, C=3, D=4, E=5

		omega_cor = this.omega_cor
		dim = this.dim

		do face = 1, 6 ! Only longitude
			do x = this.ns_xy(1), this.nf_xy(1)
				do y = this.ns_xy(2), this.nf_xy(2)
					this.f_cor(x, y, face)= 2*omega_cor*dsin(this.latlon_c(1, x, y, face)) ! function of latitude
				end do
			end do
		end do

		call this.tiles_prop(g)

		metr.r_sphere = this.r_sphere
		metr.cube_coord_c = this.cube_coord_c
		metr.latlon_c = this.latlon_c
		call metr.define(this.grid_type)

	end subroutine


end module