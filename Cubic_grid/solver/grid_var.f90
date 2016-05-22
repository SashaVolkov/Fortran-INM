module grid_var

	use sphere_geometry, Only: geometry
	use grid_generator_solver, Only: generator

implicit none

	Private
	Public :: g_var

	Type g_var

		Real(8), Allocatable :: h_dist(:, :, :)
		Real(8), Allocatable :: f_cor(:, :, :)
		Real(8), Allocatable :: points_latlon(:, :, :, :)
		Real(8), Allocatable :: square(:, :)
		Real(8), Allocatable :: triangle_area(:, :, :)
		Real(8), Allocatable :: triangle_angles(:, :, :)
		Real(8), Allocatable :: square_angles(:, :, :)
		Real(8)  omega_cor, r_sphere, g, dt, dx_min, dy_min, dx_max, dy_max, pi
		integer(4) dim, step, dim_st, rescale

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Private :: alloc => alloc
		Procedure, Private :: const_def => const_def
		Procedure, Public :: deinit => deinit

		Procedure, Private :: step_minmax => step_minmax
		Procedure, Public :: partial_c1_x => partial_c1_x
		Procedure, Public :: partial_c1_y => partial_c1_y
	End Type


CONTAINS


	subroutine init(this, geom, dim, step, omega_cor, g, dt, rescale)

		Class(g_var) :: this
		Class(geometry) :: geom
		integer(4), intent(in) :: dim, step, rescale
		real(8), intent(in) :: omega_cor, g, dt
		integer(4) x, y
		real(8) t(2), dist
		Type(generator) :: generate

		this.dim = dim;  this.step = step;  this.g = g;  this.dim_st = dim + step
		this.omega_cor = omega_cor;  this.r_sphere = geom.radius;  this.dt = dt
		this.pi = geom.pi;  this.rescale = rescale

				print '(" rad = ", f10.2, " pi = ", f10.7)', geom.radius, geom.pi

		call this.alloc()
		call generate.conformal_cubed_sphere(dim, geom.radius, rescale, this.points_latlon)
		call this.const_def(geom)
		call this.step_minmax()

	end subroutine



	subroutine alloc(this)
		Class(g_var) :: this
			Allocate(this.h_dist(1:4, -this.dim_st:this.dim_st, -this.dim_st:this.dim_st))
			Allocate(this.f_cor(-this.dim:this.dim, -this.dim:this.dim, 1:6))
			Allocate(this.points_latlon(1:2, -this.dim:this.dim, -this.dim:this.dim, 1:6))
			Allocate(this.square(-this.dim:this.dim-1, -this.dim:this.dim-1))
			Allocate(this.triangle_area(1:2, -this.dim:this.dim-1, -this.dim:this.dim-1))
			Allocate(this.triangle_angles(1:6, -this.dim:this.dim-1, -this.dim:this.dim-1))
			Allocate(this.square_angles(1:4, -this.dim:this.dim-1, -this.dim:this.dim-1))
	end subroutine



	subroutine const_def(this, g)
		Class(g_var) :: this
		Class(geometry) :: g
		real(8) dist, omega_cor, S1, S2, sphere_area
		integer(4) face, x, y, dim, step, dim_st, k
		character(8) istring


		omega_cor = this.omega_cor
		dim = this.dim;  step = this.step;  dim_st = this.dim_st


		do face = 1, 6 ! Only longitude
			do y = -dim, dim
				do x = -dim, dim
					this.f_cor(x, y, face)= 2*omega_cor*dsin(this.points_latlon(2, x, y, face))
				end do
			end do
		end do


		face = 3
		do x = -dim_st + 1, dim_st - 1
			do y = -dim_st + 1, dim_st - 1

if(y >= dim)then ! y
	call g.dist(this.points_latlon(:, 2*dim - y, x, face), this.points_latlon(:, 2*dim - y -1, x, face), dist)
	this.h_dist(1, y, x) = dist
else
	call g.dist(this.points_latlon(:, y+1, x, face), this.points_latlon(:, y, x, face), dist)
	this.h_dist(1, y, x) = dist
end if

if(x >= dim)then ! x
	call g.dist(this.points_latlon(:, y, 2*dim - x, face), this.points_latlon(:, y, 2*dim - x - 1, face), dist)
	this.h_dist(2, y, x) = dist
else
	call g.dist(this.points_latlon(:, y, x+1, face), this.points_latlon(:, y, x, face), dist)
	this.h_dist(2, y, x) = dist
end if

if(y <= -dim)then ! y
	call g.dist(this.points_latlon(:, - 2*dim - y + 1, x, face), this.points_latlon(:,  - 2*dim - y, x, face), dist)
	this.h_dist(3, y, x) = dist
else
	call g.dist(this.points_latlon(:, y, x, face), this.points_latlon(:, y-1, x, face), dist)
	this.h_dist(3, y, x) = dist
end if

if(x <= -dim)then ! x
	call g.dist(this.points_latlon(:, y, - 2*dim - x + 1, face), this.points_latlon(:, y,  - 2*dim - x, face), dist)
	this.h_dist(4, y, x) = dist
else
	call g.dist(this.points_latlon(:, y, x, face), this.points_latlon(:, y, x-1, face), dist)
	this.h_dist(4, y, x) = dist
end if

			end do
		end do

		sphere_area = 0



		if (this.rescale == 1) then
			istring = '_tan'
		else if (this.rescale == 0) then
			istring = '_simple'
		else if (this.rescale == 2) then
			istring = '_exp'
		end if

		open (20, file = '/home/sasha/Fortran/Cubic_grid/solver/datFiles/angle'//trim(istring)//'.dat')

		do x = -dim, dim-1
			do y = -dim, dim-1

				call g.triangle(this.points_latlon(:, y, x, face), this.points_latlon(:, y, x+1, face), this.points_latlon(:, y+1, x, face), S1, this.triangle_angles(1:3, y, x))
				call g.triangle(this.points_latlon(:, y+1, x+1, face), this.points_latlon(:, y, x+1, face), this.points_latlon(:, y+1, x, face), S2, this.triangle_angles(4:6, y, x))

				this.square(x, y) = S1 + S2
				this.triangle_area(1, x, y) = S1
				this.triangle_area(2, x, y) = S2

				sphere_area = sphere_area + this.square(x,y)

				this.square_angles(1, x, y) = this.triangle_angles(2, y, x)
				this.square_angles(2, x, y) = this.triangle_angles(3, y, x) + this.triangle_angles(6, y, x)
				this.square_angles(3, x, y) = this.triangle_angles(5, y, x)
				this.square_angles(4, x, y) = this.triangle_angles(1, y, x) + this.triangle_angles(4, y, x)

				do k=1,4
					write(20,*) this.square_angles(k, x, y)*180d0/this.pi
				end do

			end do
		end do

		close(20)

		print '(" sphere_area = ", f20.2)', sphere_area*6d0

	end subroutine



	subroutine deinit(this)
		Class(g_var) :: this
			if (Allocated(this.h_dist)) Deallocate(this.h_dist)
			if (Allocated(this.f_cor)) Deallocate(this.f_cor)
			if (Allocated(this.points_latlon)) Deallocate(this.points_latlon)
			if (Allocated(this.square)) Deallocate(this.square)
			if (Allocated(this.triangle_area)) Deallocate(this.triangle_area)
			if (Allocated(this.triangle_angles)) Deallocate(this.triangle_angles)
			if (Allocated(this.square_angles)) Deallocate(this.square_angles)
	end subroutine



	real(8) function partial_c1_x(this, fun, x, y)
		Class(g_var) :: this
		real(8), intent(in) :: fun(1:3)
		integer, intent(in) :: x, y
		partial_c1_x = (fun(3) - fun(1))/(this.h_dist(4, x, y) + this.h_dist(2, x, y))

	end function



	real(8) function partial_c1_y(this, fun, x, y)
		Class(g_var) :: this
		real(8), intent(in) :: fun(1:3)
		integer, intent(in) :: x, y
		partial_c1_y = (fun(3) - fun(1))/(this.h_dist(3, x, y) + this.h_dist(1, x, y))

	end function



subroutine step_minmax(this)
	Class(g_var) :: this
	real(8) dim

	dim = this.dim
	this.dx_min = MINVAL(this.h_dist(2,-dim:dim,-dim:dim))
	this.dy_min = MINVAL(this.h_dist(1,-dim:dim,-dim:dim))
	this.dx_max = MAXVAL(this.h_dist(2,-dim:dim,-dim:dim))
	this.dy_max = MAXVAL(this.h_dist(1,-dim:dim,-dim:dim))

end subroutine




end module