		! ____________________
		! | 1, 2d     2d, 2d |
		! |                  |
		! |                  |     ___
		! |       face       |   __|6|____
		! |    orientation   |   |5|2|3|4| - faces
		! |                  |     |1|
		! |                  |
		! | 1,1        2d, 1 |
		! --------------------

module grid_var

	use sphere_geometry, Only: geometry
	use grid_generator_solver, Only: generator
	use parallel_cubic, Only: parallel
	use metrics, Only: metric

implicit none

	Private
	Public :: g_var

	Type g_var

		Real(8), Allocatable :: f_cor(:, :, :)
		Real(8), Allocatable :: latlon_c(:, :, :, :)
		Real(8), Allocatable :: cube_coord_c(:, :, :)
		Real(8), Allocatable :: latlon(:, :, :, :)

		Real(8), Allocatable :: square(:, :)
		Real(8), Allocatable :: real_dist(:, :)
		Real(8), Allocatable :: triangle_area(:, :, :)
		Real(8), Allocatable :: triangle_angles(:, :, :)
		Real(8), Allocatable :: square_angles(:, :, :)

		Real(8) :: four_order_const(5)   ! "Compact finite difference schemes on non-uniform meshes" Gamet et al. 1999 
		Real(8) ::  omega_cor, r_sphere, g, dt, dx_min, dy_min, dx_max, dy_max, pi, delta_on_cube
		integer(4) dim, step, rescale, ns_xy(2), nf_xy(2), grid_type
		integer(4) first_x, first_y, last_x, last_y, snd_xy(6, 4, 2), rcv_xy(6, 4, 2)

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Private :: alloc => alloc
		Procedure, Public :: deinit => deinit
		Procedure, Private :: const_def => const_def
		Procedure, Private :: tiles_prop => tiles_prop
	End Type


CONTAINS


	subroutine init(this, geom, paral, metr, omega_cor, g, dt, rescale,grid_type)

		Class(g_var) :: this
		Class(geometry) :: geom
		Class(parallel) :: paral
		Class(metric) :: metr
		integer(4), intent(in) :: rescale, grid_type
		real(8), intent(in) :: omega_cor, g, dt
		integer(4) x, y
		real(8) t(2), dist
		Type(generator) :: generate

		this.dim = paral.dim;  this.step = paral.step;  this.g = g
		this.omega_cor = omega_cor;  this.r_sphere = geom.radius;  this.dt = dt
		this.pi = geom.pi;  this.rescale = rescale;  this.grid_type = grid_type

		this.ns_xy(:) = paral.ns_xy(:);  this.nf_xy(:) = paral.nf_xy(:);
		this.first_x = paral.first_x;  this.first_y = paral.first_y
		this.last_x = paral.last_x;  this.last_y = paral.last_y

		this.snd_xy(:,:,:) = paral.snd_xy(:,:,:)
		this.rcv_xy(:,:,:) = paral.rcv_xy(:,:,:)


		call this.alloc()
		if(grid_type == 0) then
			call generate.conformal_cubed_sphere(this.dim, this.step, this.r_sphere, rescale, this.cube_coord_c, this.latlon_c, this.latlon)
		else if(grid_type == 1)then
			call generate.equiangular_cubed_sphere(this.dim, this.step, this.cube_coord_c, this.latlon_c, this.latlon)
		end if
		call this.const_def(geom, metr)

	end subroutine



	subroutine alloc(this)
		Class(g_var) :: this
		integer(4) f_x, f_y, l_x, l_y, dim, step, f, l

		f_x = this.first_x;  l_x = this.last_x;  f_y = this.first_y;  l_y = this.last_y
		dim = this.dim;  step = this.step;  f = 1-step; l = 2*dim + step

		Allocate(this.f_cor(f_x:l_x , f_y:l_y, 1:6))
		Allocate(this.latlon_c(1:2, f:l, f:l, 1:6))
		Allocate(this.cube_coord_c(1:2, f:l , f:l))
		Allocate(this.latlon(1:2, f:l+1 , f:l+1, 1:6))

		Allocate(this.square(1:2*dim, 1:2*dim))
		Allocate(this.real_dist(this.ns_xy(1):this.nf_xy(1), this.ns_xy(2):this.nf_xy(2)))
		Allocate(this.triangle_area(1:2, 1:2*dim, 1:2*dim))
		Allocate(this.triangle_angles(1:6, 1:2*dim, 1:2*dim))
		Allocate(this.square_angles(1:4, 1:2*dim, 1:2*dim))

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



	subroutine tiles_prop(this, g)
		Class(g_var) :: this
		Class(geometry) :: g
		real(8) dist, omega_cor, S1, S2, sphere_area
		integer(4) x, y, dim, k, step
		character(8) istring


		omega_cor = this.omega_cor
		dim = this.dim;  step = this.step
		sphere_area = 0

		this.delta_on_cube = (this.cube_coord_c(1, dim, dim) - this.cube_coord_c(1, dim-1, dim))*this.r_sphere

		do x = 1, 2*dim
			do y = 1, 2*dim

	call g.triangle(this.latlon(:, y, x, 2), this.latlon(:, y, x+1, 2), this.latlon(:, y+1, x, 2), S1, this.triangle_angles(1:3, y, x))
	call g.triangle(this.latlon(:, y+1, x+1, 2), this.latlon(:, y, x+1, 2), this.latlon(:, y+1, x, 2), S2, this.triangle_angles(4:6, y, x))

	this.square(x, y) = S1 + S2
	this.triangle_area(1, x, y) = S1
	this.triangle_area(2, x, y) = S2

	sphere_area = sphere_area + this.square(x,y)

	this.square_angles(1, x, y) = this.triangle_angles(2, y, x)
	this.square_angles(2, x, y) = this.triangle_angles(3, y, x) + this.triangle_angles(6, y, x)
	this.square_angles(3, x, y) = this.triangle_angles(5, y, x)
	this.square_angles(4, x, y) = this.triangle_angles(1, y, x) + this.triangle_angles(4, y, x)

			end do
		end do

		do x = this.ns_xy(1), this.nf_xy(1)
			do y = this.ns_xy(2), this.nf_xy(2)
				this.real_dist(x,y) = g.dist(this.latlon_c(:, y, x, 2), this.latlon_c(:, y, x-1, 2))
			end do
		end do

	end subroutine



end module