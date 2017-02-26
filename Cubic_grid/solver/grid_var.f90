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
	use omp_lib

implicit none

	Private
	Public :: g_var

	Type g_var

		Real(8), Allocatable :: f_cor(:, :, :)
		Real(8), Allocatable :: latlon_c(:, :, :, :)
		Real(8), Allocatable :: cube_coord_c(:, :, :)
		Real(8), Allocatable :: latlon(:, :, :, :)

		Real(8), Allocatable :: square(:, :)
		Real(8), Allocatable :: Real_dist(:, :)
		Real(8), Allocatable :: triangle_area(:, :, :)
		Real(8), Allocatable :: triangle_angles(:, :, :)
		Real(8), Allocatable :: square_angles(:, :, :)

		Real(8) :: four_order_const(5)   ! "Compact finite difference schemes on non-uniform meshes" Gamet et al. 1999 
		Real(8) ::  omega_cor, r_sphere, dx_min, dy_min, dx_max, dy_max, pi, delta_on_cube, max_to_min, min
		Integer(4) dim, step, rescale, ns_xy(2), nf_xy(2), grid_type, Neighbours_face(6, 4), id
		Integer(4) first_x, first_y, last_x, last_y, snd_xy(6, 4, 2), rcv_xy(6, 4, 2)

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Private :: alloc => alloc
		Procedure, Public :: deinit => deinit
		Procedure, Private :: tiles_prop => tiles_prop
	End Type


CONTAINS


	Subroutine init(this, geom, paral, rescale,grid_type)

		Class(g_var) :: this
		Class(geometry) :: geom
		Class(parallel) :: paral
		Integer(4), intent(in) :: rescale, grid_type
		Integer(4) x, y
		Real(8) t(2), dist
		Type(generator) :: generate

		this.dim = paral.dim;  this.step = paral.step;  this.id = paral.id
		this.r_sphere = geom.radius
		this.pi = geom.pi;  this.rescale = rescale;  this.grid_type = grid_type

		this.ns_xy(:) = paral.ns_xy(:);  this.nf_xy(:) = paral.nf_xy(:);
		this.first_x = paral.first_x;  this.first_y = paral.first_y
		this.last_x = paral.last_x;  this.last_y = paral.last_y

		this.snd_xy = paral.snd_xy;  this.rcv_xy = paral.rcv_xy
		this.Neighbours_face = paral.Neighbours_face


		call this.alloc()
		if(grid_type == 0) then
			call generate.conformal_cubed_sphere(this.dim, this.step, this.r_sphere, rescale, this.cube_coord_c, this.latlon_c, this.latlon)
		else if(grid_type == 1)then
			call generate.equiangular_cubed_sphere(this.dim, this.step, this.cube_coord_c, this.latlon_c, this.latlon)
		end if
		call this.tiles_prop(geom)

	end Subroutine



	Subroutine alloc(this)
		Class(g_var) :: this
		Integer(4) f_x, f_y, l_x, l_y, dim, step, f, l

		f_x = this.first_x;  l_x = this.last_x;  f_y = this.first_y;  l_y = this.last_y
		dim = this.dim;  step = this.step;  f = 1-2*step; l = 2*dim + 2*step

		Allocate(this.f_cor(f_x:l_x , f_y:l_y, 6))
		Allocate(this.latlon_c(2, f:l, f:l, 6))
		Allocate(this.cube_coord_c(2, f:l , f:l))
		Allocate(this.latlon(2, f:l+1 , f:l+1, 6))

		Allocate(this.square(1:2*dim, 1:2*dim))
		Allocate(this.Real_dist(1:2*dim, 1:2*dim))
		Allocate(this.triangle_area(1:2, 1:2*dim, 1:2*dim))
		Allocate(this.triangle_angles(1:6, 1:2*dim, 1:2*dim))
		Allocate(this.square_angles(1:4, 1:2*dim, 1:2*dim))

	end Subroutine



	Subroutine deinit(this)
		Class(g_var) :: this
		if (Allocated(this.latlon_c)) Deallocate(this.latlon_c)
		if (Allocated(this.cube_coord_c)) Deallocate(this.cube_coord_c)
		if (Allocated(this.latlon)) Deallocate(this.latlon)

		if (Allocated(this.square)) Deallocate(this.square)
		if (Allocated(this.triangle_area)) Deallocate(this.triangle_area)
		if (Allocated(this.triangle_angles)) Deallocate(this.triangle_angles)
		if (Allocated(this.square_angles)) Deallocate(this.square_angles)
	end Subroutine





	Subroutine tiles_prop(this, g)
		Class(g_var) :: this
		Class(geometry) :: g
		Real(8) dist, omega_cor, S1, S2, sphere_area
		Integer(4) x, y, dim, k, step, id, ier
		character(8) istring


		omega_cor = this.omega_cor
		dim = this.dim;  step = this.step
		sphere_area = 0


		if(this.grid_type == 0) then
			this.delta_on_cube = 1d0/Real(dim, 8)*this.r_sphere
		else if(this.grid_type == 1)then
			this.delta_on_cube = (this.cube_coord_c(1, dim, dim) - this.cube_coord_c(1, dim-1, dim))*this.r_sphere
		end if

		!$OMP PARALLEL PRIVATE(y, x, sphere_area)
		!$OMP DO

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

		!$OMP END DO
		!$OMP DO

		do x = 2, 2*dim
			do y = 1, 2*dim
				this.Real_dist(x,y) = g.dist(this.latlon_c(:, y, x, 2), this.latlon_c(:, y, x-1, 2))
			end do
		end do

		!$OMP END DO
		!$OMP END PARALLEL

		this.min = minval(this.Real_dist(2*dim-2:2*dim, 2*dim-2:2*dim))
		this.max_to_min = maxval(this.Real_dist(2*dim-2:2*dim, 2*dim-2:2*dim))/this.min

	end Subroutine



end module