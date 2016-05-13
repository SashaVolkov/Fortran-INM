module grid_var

	use geometry
	use grid_generator_solver, Only: generator

implicit none

	Private
	Public :: g_var

	Type g_var

		Real(8), Allocatable :: h_dist(:, :, :)
		Real(8), Allocatable :: f_cor(:, :, :)
		Real(8), Allocatable :: grid_points_latlon(:, :, :, :)
! 		Real(8), Allocatable :: grid_points_xy(:, :, :, :)
! 		Real(8), Allocatable :: grid_points_xy_tan(:, :, :, :)
! 		Real(8), Allocatable :: grid_points_xyz(:, :, :, :)
		Real(8)  omega_cor, r_sphere, g, dt, dx_min, dy_min, dx_max, dy_max, pi
		integer(4) dim, step, dim_st

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



	subroutine init(this, dim, step, omega_cor, r_sphere, g, dt)

		Class(g_var) :: this
		integer(4), intent(in) :: dim, step ! dimension
		real(8), intent(in) :: omega_cor, r_sphere, g, dt
		integer(4) face_idx, x, y
		real(8) t(2), dist
		Type(generator) :: generate

		this.dim = dim;  this.step = step;  this.g = g;  this.dim_st = dim + step
		this.omega_cor = omega_cor;  this.r_sphere = r_sphere;  this.dt = dt


		this.pi = 314159265358979323846d-20

		print '(" allocate")'
		call this.alloc()
		print '(" generate")'
		call generate.conformal_cubed_sphere(dim, r_sphere, this.grid_points_latlon)
		print '(" define")'
		call this.const_def()
		call this.step_minmax()

	end subroutine



	subroutine alloc(this)

		Class(g_var) :: this

		Allocate(this.h_dist(1:4, -this.dim_st:this.dim_st, -this.dim_st:this.dim_st))
		Allocate(this.f_cor(-this.dim:this.dim, -this.dim:this.dim, 1:6))
		Allocate(this.grid_points_latlon(1:2, -this.dim:this.dim, -this.dim:this.dim, 1:6))
! 		Allocate(this.grid_points_xy(1:2, -this.dim_st:this.dim_st, -this.dim_st:this.dim_st, 1:6))
! 		Allocate(this.grid_points_xy_tan(1:2, -this.dim_st:this.dim_st, -this.dim_st:this.dim_st, 1:6))

	end subroutine



	subroutine const_def(this)
		Class(g_var) :: this
		! integer(4), intent(in) :: dim, step ! dimension
		real(8) grid_points(1:2, -this.dim:this.dim, -this.dim:this.dim)
		real(8) dist, omega_cor, r_sphere
		integer(4) face_idx, x, y, dim, step, dim_st

		omega_cor = this.omega_cor;  r_sphere = this.r_sphere
		dim = this.dim;  step = this.step;  dim_st = this.dim_st
		grid_points = this.grid_points_latlon(:, :, :, 2)

		do face_idx = 1, 6 ! Only longitude
			do y = -dim, dim
				do x = -dim, dim
					this.f_cor(x, y, face_idx)= 2*omega_cor*dsin(this.grid_points_latlon(2, x, y, face_idx))
				end do
			end do
		end do


		face_idx = 2
		do x = -dim, dim
			do y = -dim, dim

if(y+1 > dim)then ! y
	call distance_sphere(r_sphere, grid_points(:, y, x), grid_points(:, y-1, x), dist)
	this.h_dist(1, y, x) = dist
else
	call distance_sphere(r_sphere, grid_points(:, y+1, x), grid_points(:, y, x), dist)
	this.h_dist(1, y, x) = dist
end if

if(x+1 > dim)then ! x
	call distance_sphere(r_sphere, grid_points(:, y, x), grid_points(:, y, x-1), dist)
	this.h_dist(2, y, x) = dist
else
	call distance_sphere(r_sphere, grid_points(:, y, x+1), grid_points(:, y, x), dist)
	this.h_dist(2, y, x) = dist
end if

if(y-1 < -dim)then ! y
	call distance_sphere(r_sphere, grid_points(:, y+1, x), grid_points(:, y, x), dist)
	this.h_dist(3, y, x) = dist
else
	call distance_sphere(r_sphere, grid_points(:, y, x), grid_points(:, y-1, x), dist)
	this.h_dist(3, y, x) = dist
end if

if(x-1 < -dim)then ! x
	call distance_sphere(r_sphere, grid_points(:, y, x+1), grid_points(:, y, x), dist)
	this.h_dist(4, y, x) = dist
else
	call distance_sphere(r_sphere, grid_points(:, y, x), grid_points(:, y, x-1), dist)
	this.h_dist(4, y, x) = dist
end if

			end do
		end do

	end subroutine



	subroutine deinit(this)
		Class(g_var) :: this

		if (Allocated(this.h_dist)) Deallocate(this.h_dist)
		if (Allocated(this.f_cor)) Deallocate(this.f_cor)
		if (Allocated(this.grid_points_latlon)) Deallocate(this.grid_points_latlon)
! 		if (Allocated(this.grid_points_xy)) Deallocate(this.grid_points_xy)
! 		if (Allocated(this.grid_points_xy_tan)) Deallocate(this.grid_points_xy_tan)

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