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

implicit none

	Private
	Public :: g_var

	Type g_var

		Real(8), Allocatable :: x_dist(:, :)
		Real(8), Allocatable :: y_dist(:, :)
		Real(8), Allocatable :: G_sqr(:, :)
		Real(8), Allocatable :: G_tensor(:, :, :, :)
		Real(8), Allocatable :: G_inverse(:, :, :, :)
		Real(8), Allocatable :: rho(:, :)
		Real(8), Allocatable :: To_sph_coord(:, :, :, :, :)
		Real(8), Allocatable :: From_sph_coord(:, :, :, :, :)

		Real(8), Allocatable :: f_cor(:, :, :)
		Real(8), Allocatable :: latlon_c(:, :, :, :)
		Real(8), Allocatable :: equiang_c(:, :, :, :)
		Real(8), Allocatable :: latlon(:, :, :, :)

		Real(8), Allocatable :: square(:, :)
		Real(8), Allocatable :: triangle_area(:, :, :)
		Real(8), Allocatable :: triangle_angles(:, :, :)
		Real(8), Allocatable :: square_angles(:, :, :)

		Real(8), Allocatable :: four_order_const_x(:, :, :)   ! "Compact finite difference schemes on non-uniform meshes" Gamet et al. 1999 
		Real(8), Allocatable :: four_order_const_y(:, :, :)
		Real(8) ::  omega_cor, r_sphere, g, dt, dx_min, dy_min, dx_max, dy_max, pi
		integer(4) dim, step, rescale, ns_xy(2), nf_xy(2), grid_type
		integer(4) first_x, first_y, last_x, last_y

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Private :: alloc => alloc
		Procedure, Public :: deinit => deinit
		Procedure, Private :: const_def => const_def
		Procedure, Private :: step_minmax => step_minmax

		Procedure, Private :: transformation_matrix_equiang => transformation_matrix_equiang
		Procedure, Private :: transformation_sph_equiang => transformation_sph_equiang
		Procedure, Private :: tiles_prop => tiles_prop
		Procedure, Private :: derivat_4order => derivat_4order
	End Type


CONTAINS


	subroutine init(this, geom, paral, omega_cor, g, dt, rescale,grid_type)

		Class(g_var) :: this
		Class(geometry) :: geom
		Class(parallel) :: paral
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


				! print '(" rad = ", f10.2, " pi = ", f10.7)', geom.radius, geom.pi

		call this.alloc()
		if(grid_type == 0) then
			call generate.conformal_cubed_sphere(this.dim, this.r_sphere, rescale, this.latlon_c, this.latlon)
		else if(grid_type == 1)then
			call generate.equiangular_cubed_sphere(this.dim, this.step, this.equiang_c, this.latlon_c, this.latlon)
		end if
		call this.const_def(geom)
		call this.step_minmax()

	end subroutine



	subroutine alloc(this)
		Class(g_var) :: this
		integer(4) f_x, f_y, l_x, l_y, dim, step, f, l

		f_x = this.first_x;  l_x = this.last_x;  f_y = this.first_y;  l_y = this.last_y
		dim = this.dim;  step = this.step;  f = 1-step; l = 2*dim + step

			Allocate(this.x_dist(f:l, f:l))
			Allocate(this.y_dist(f:l, f:l))
			Allocate(this.G_sqr(f_x:l_x , f_y:l_y))
			Allocate(this.G_tensor(f_x:l_x , f_y:l_y, 1:2, 1:2))
			Allocate(this.G_inverse(f_x:l_x , f_y:l_y, 1:2, 1:2))
			Allocate(this.rho(f_x:l_x , f_y:l_y))
			Allocate(this.To_sph_coord(2, 2, f_x:l_x , f_y:l_y, 6))
			Allocate(this.From_sph_coord(2, 2, f_x:l_x , f_y:l_y, 6))

			Allocate(this.f_cor(f_x:l_x , f_y:l_y, 1:6))
			Allocate(this.latlon_c(1:2, f:l, f:l, 1:6))
			Allocate(this.equiang_c(1:2, f:l , f:l, 1:6))
			Allocate(this.latlon(1:2, f:l+1 , f:l+1, 1:6))

			Allocate(this.square(1:2*dim, 1:2*dim))
			Allocate(this.triangle_area(1:2, 1:2*dim, 1:2*dim))
			Allocate(this.triangle_angles(1:6, 1:2*dim, 1:2*dim))
			Allocate(this.square_angles(1:4, 1:2*dim, 1:2*dim))

			Allocate(this.four_order_const_x(1:5, f_x:l_x , f_y:l_y))
			Allocate(this.four_order_const_y(1:5, f_x:l_x , f_y:l_y))
	end subroutine



	subroutine deinit(this)
		Class(g_var) :: this
			if (Allocated(this.x_dist)) Deallocate(this.x_dist)
			if (Allocated(this.y_dist)) Deallocate(this.y_dist)
			if (Allocated(this.G_sqr)) Deallocate(this.G_sqr)
			if (Allocated(this.G_tensor)) Deallocate(this.G_tensor)
			if (Allocated(this.G_inverse)) Deallocate(this.G_inverse)
			if (Allocated(this.rho)) Deallocate(this.rho)
			if (Allocated(this.To_sph_coord)) Deallocate(this.To_sph_coord)
			if (Allocated(this.From_sph_coord)) Deallocate(this.From_sph_coord)

			if (Allocated(this.f_cor)) Deallocate(this.f_cor)
			if (Allocated(this.latlon_c)) Deallocate(this.latlon_c)
			if (Allocated(this.equiang_c)) Deallocate(this.equiang_c)
			if (Allocated(this.latlon)) Deallocate(this.latlon)

			if (Allocated(this.square)) Deallocate(this.square)
			if (Allocated(this.triangle_area)) Deallocate(this.triangle_area)
			if (Allocated(this.triangle_angles)) Deallocate(this.triangle_angles)
			if (Allocated(this.square_angles)) Deallocate(this.square_angles)

			if (Allocated(this.four_order_const_x)) Deallocate(this.four_order_const_x)
			if (Allocated(this.four_order_const_y)) Deallocate(this.four_order_const_y)
	end subroutine



	subroutine const_def(this, g)
		Class(g_var) :: this
		Class(geometry) :: g
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
		call this.transformation_matrix_equiang()
		if(grid_type == 0) then
			call this.transformation_matrix_conf()
		else if(grid_type == 1)then
			call this.transformation_matrix_equiang()
		end if
		if(this.step > 1) call this.derivat_4order()

	end subroutine




	subroutine derivat_4order(this)
		Class(g_var) :: this
		real(8) dist, sphere_area, h(-1:2)
		integer(4) face, x, y, dim, step
		integer(4), parameter :: A =1, B=2, C=3, D=4, E=5

		dim = this.dim;  step = this.step

		do x = this.ns_xy(1), this.nf_xy(1)
			do y = this.ns_xy(2), this.nf_xy(2) ! Gamet et al. 1999 Apendix A. Approx. of derivat.

				h(-1) = this.x_dist(x-1, y);  h(0) = this.x_dist(x, y)
				h(1) = this.x_dist(x+1, y);  h(2) = this.x_dist(x+2, y)

this.four_order_const_x( A, x, y) = ( h(1) + h(2) )*( h(-1)*h(0) + h(0)**2 )/( h(1)*h(2)*( h(0) + h(1) )*( h(-1) + h(0) + h(1) ) )

this.four_order_const_x( B, x, y) = - ( h(-1) + h(0) )*( h(1)*h(2) + h(1)**2 )/( h(-1)*h(0)*( h(0) + h(1) )*( h(0) + h(1) + h(2) ) )

this.four_order_const_x( C, x, y) = - ( h(-1) + h(0) )*h(0)*h(1)/( h(2)*( h(1) + h(2) )*( h(0) + h(1) + h(2) )*( h(-1) + h(0) + h(1) + h(2) ) )

this.four_order_const_x( D, x, y) = ( h(1) + h(2) )*h(0)*h(1)/( h(-1)*( h(-1) + h(0) )*( h(-1) + h(0) + h(1) )*( h(-1) + h(0) + h(1) + h(2) ) )

this.four_order_const_x( E, x, y) = - ( this.four_order_const_x( A, x, y) + this.four_order_const_x( B, x, y) + this.four_order_const_x( C, x, y) + this.four_order_const_x( D, x, y) )

				end do
			end do

			do x = this.ns_xy(1), this.nf_xy(1)
				do y = this.ns_xy(2), this.nf_xy(2)

				h(-1) = this.y_dist(x, y-1);  h(0) = this.y_dist(x, y)
				h(1) = this.y_dist(x, y+1);  h(2) = this.y_dist(x, y+2)

this.four_order_const_y( A, x, y) = ( h(1) + h(2) )*( h(-1)*h(0) + h(0)**2 )/( h(1)*h(2)*( h(0) + h(1) )*( h(-1) + h(0) + h(1) ) )

this.four_order_const_y( B, x, y) = - ( h(-1) + h(0) )*( h(1)*h(2) + h(1)**2 )/( h(-1)*h(0)*( h(0) + h(1) )*( h(0) + h(1) + h(2) ) )

this.four_order_const_y( C, x, y) = - ( h(-1) + h(0) )*h(0)*h(1)/( h(2)*( h(1) + h(2) )*( h(0) + h(1) + h(2) )*( h(-1) + h(0) + h(1) + h(2) ) )

this.four_order_const_y( D, x, y) = ( h(1) + h(2) )*h(0)*h(1)/( h(-1)*( h(-1) + h(0) )*( h(-1) + h(0) + h(1) )*( h(-1) + h(0) + h(1) + h(2) ) )

this.four_order_const_y( E, x, y) = - ( this.four_order_const_y( A, x, y) + this.four_order_const_y( B, x, y) + this.four_order_const_y( C, x, y) + this.four_order_const_y( D, x, y) )

				end do
			end do


	end subroutine



	subroutine transformation_matrix_equiang(this)
		Class(g_var) :: this
		real(8) x_1, x_2, g_coef, g_inv_coef
		integer(4) x, y, face


		do y = this.first_y, this.last_y
			do x = this.first_x, this.last_x

			x_1 = this.equiang_c(1, x, y, 2)
			x_2 = this.equiang_c(2, x, y, 2)

			this.rho(x, y) = dsqrt(1 + (dtan(x_1))**2 + (dtan(x_2))**2)
			this.G_sqr(x, y) = 1.0/(this.rho(x, y)**3 * ((dcos(x_1))**2) * (dcos(x_2)**2))
			g_coef = 1.0/((this.rho(x, y)**4) * ((dcos(x_1))**2) * (dcos(x_2)**2))
			g_inv_coef = ((this.rho(x, y)**2) * ((dcos(x_1))**2) * (dcos(x_2)**2))

			this.G_tensor(x, y, 1, 1) = g_coef * (1 + (dtan(x_1))**2)
			this.G_tensor(x, y, 1, 2) = g_coef * (- (dtan(x_1))*(dtan(x_2)))
			this.G_tensor(x, y, 2, 1) = g_coef * (- (dtan(x_1))*(dtan(x_2)))
			this.G_tensor(x, y, 2, 2) = g_coef * (1 + (dtan(x_2))**2)

			this.G_inverse(x, y, 1, 1) = g_inv_coef * (1 + (dtan(x_2))**2)
			this.G_inverse(x, y, 1, 2) = g_inv_coef * (dtan(x_1))*(dtan(x_2))
			this.G_inverse(x, y, 2, 1) = g_inv_coef * (dtan(x_1))*(dtan(x_2))
			this.G_inverse(x, y, 2, 2) = g_inv_coef * (1 + (dtan(x_1))**2)

			end do
		end do

		call this.transformation_sph_equiang()

	end subroutine



	subroutine transformation_matrix_conf(this)
		Class(g_var) :: this
		integer(4) x, y, face

		do y = this.first_y, this.last_y
			do x = this.first_x, this.last_x

			this.G_sqr(x, y) = 1.0

			this.G_tensor(x, y, 1, 1) = 1.0
			this.G_tensor(x, y, 1, 2) = 0.0
			this.G_tensor(x, y, 2, 1) = 0.0
			this.G_tensor(x, y, 2, 2) = 1.0

			this.G_inverse(x, y, 1, 1) = 1.0
			this.G_inverse(x, y, 1, 2) = 0.0
			this.G_inverse(x, y, 2, 1) = 0.0
			this.G_inverse(x, y, 2, 2) = 1.0

			end do
		end do

	end subroutine



	subroutine transformation_sph_equiang(this)
		Class(g_var) :: this
		real(8) x_1, x_2, g_coef, s(6)
		integer(4) x, y, face, delta
		s(1) = - 1d0;  s(6) = 1d0


		do face = 2,5
			do y = this.first_y, this.last_y
				do x = this.first_x, this.last_x
					delta = this.rho(x,y)
					x_1 = this.equiang_c(1, x, y, face)
					x_2 = this.equiang_c(2, x, y, face)

					this.From_sph_coord(1,1,x,y,face) = 1d0
					this.From_sph_coord(1,2,x,y,face) = 0d0
					this.From_sph_coord(2,1,x,y,face) = dtan(x_1)*dtan(x_2)*((dcos(x_2))**2)
					this.From_sph_coord(2,2,x,y,face) = (delta**2)*((dcos(x_2))**2)*abs(dcos(x_1))

					this.To_sph_coord(1,1,x,y,face) = 1d0
					this.To_sph_coord(1,2,x,y,face) = 0d0
					this.To_sph_coord(2,1,x,y,face) = dtan(x_1)*dtan(x_2)/(abs(dcos(x_1))*(delta**2))
					this.To_sph_coord(2,2,x,y,face) = 1d0/this.From_sph_coord(2,2,x,y,face)
				end do
			end do
		end do

		do face = 1, 6, 5
			do y = this.first_y, this.last_y
				do x = this.first_x, this.last_x
					delta = this.rho(x,y)
					x_1 = this.equiang_c(1, x, y, face)
					x_2 = this.equiang_c(2, x, y, face)

					this.From_sph_coord(1,1,x,y,face) = -s(face)*dtan(x_2)*((dcos(x_1))**2)
					this.From_sph_coord(1,2,x,y,face) = -s(face)*(delta**2)*dtan(x_1)*((dcos(x_1))**2)/dsqrt(dtan(x_1)**2 + dtan(x_2)**2)
					this.From_sph_coord(2,1,x,y,face) = s(face)*dtan(x_1)*((dcos(x_2))**2)
					this.From_sph_coord(2,2,x,y,face) = -s(face)*(delta**2)*dtan(x_2)*((dcos(x_2))**2)/dsqrt(dtan(x_1)**2 + dtan(x_2)**2)

					this.To_sph_coord(1,1,x,y,face) = - s(face)*dtan(x_2)/(((dcos(x_1))**2) * (dtan(x_1)**2 + dtan(x_2)**2))
					this.To_sph_coord(1,2,x,y,face) = s(face)*dtan(x_1)/(((dcos(x_2))**2) * (dtan(x_1)**2 + dtan(x_2)**2))
					this.To_sph_coord(2,1,x,y,face) = - s(face)*dtan(x_1)/(((dcos(x_1))**2) * dsqrt(dtan(x_1)**2 + dtan(x_2)**2) * (delta**2))
					this.To_sph_coord(2,2,x,y,face) = - s(face)*dtan(x_2)/(((dcos(x_2))**2) * dsqrt(dtan(x_1)**2 + dtan(x_2)**2) * (delta**2))
				end do
			end do
		end do

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

		do y = 1-step, 2*dim + step
			do x = 2-step, 2*dim + step
				if(this.grid_type == 0) then

					this.x_dist(x, y) = g.dist(this.latlon_c(:, x, y, 2), this.latlon_c(:, x-1, y, 2))*this.r_sphere
					this.y_dist(y, x) = this.x_dist(x, y)

				else if(this.grid_type == 1)then

					this.x_dist(x, y) = (this.equiang_c(1, x, y, 2) - this.equiang_c(1, x-1, y, 2))*this.r_sphere
					this.y_dist(y, x) = this.x_dist(x, y)

				end if
			end do
		end do


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


	end subroutine



subroutine step_minmax(this)
	Class(g_var) :: this
	real(8) dim

	dim = this.dim
	this.dx_min = MINVAL(this.x_dist(2:2*dim, 1:2*dim))
	this.dy_min = MINVAL(this.y_dist(1:2*dim, 2:2*dim))
	this.dx_max = MAXVAL(this.x_dist(2:2*dim, 1:2*dim))
	this.dy_max = MAXVAL(this.y_dist(1:2*dim, 2:2*dim))

end subroutine



end module