module grid_var

use geometry

implicit none

	Private
	Public :: g_var

	Type g_var

		Real(8), Allocatable :: h_dist(:, :, :)
		Real(8), Allocatable :: f_cor(:, :, :)
		real(8) g
		integer(4) dim, step, dim_st

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Private :: alloc => alloc
		Procedure, Private :: const_def => const_def
		Procedure, Public :: deinit => deinit
		Procedure, Public :: partial_c1_x => partial_c1_x
		Procedure, Public :: partial_c1_y => partial_c1_y
	End Type


CONTAINS



	subroutine init(this, grid_points, dim, step, omega_cor, r_sphere, g)

		Class(g_var) :: this
		integer(4), intent(in) :: dim, step ! dimension
		real(8), intent(in) :: grid_points(1:2, -dim:dim, -dim:dim, 1:6), omega_cor, r_sphere, g

		this.dim = dim;  this.step = step;  this.g = g;
		this.dim_st = dim + step

		call this.alloc()
		call this.const_def(grid_points, dim, step, omega_cor, r_sphere)

	! print '("  Grid step real = ", f10.3, " m")', this.h_dist(3, 300, 200, 2)

	end subroutine



	subroutine alloc(this)

		Class(g_var) :: this

		Allocate(this.h_dist(4*this.step, -this.dim_st:this.dim_st, -this.dim_st:this.dim_st))
		Allocate(this.f_cor(-this.dim_st:this.dim_st, -this.dim_st:this.dim_st, 1:6))

	end subroutine



	subroutine const_def(this, grid_points, dim, step, omega_cor, r_sphere)
		Class(g_var) :: this
		integer(4), intent(in) :: dim, step ! dimension
		real(8), intent(in) :: grid_points(1:2, -dim:dim, -dim:dim, 1:6), omega_cor, r_sphere
		real(8) dist
		integer(4) face_idx, x, y


		do face_idx = 1, 6 ! Only longitude
			do y = -this.dim_st, this.dim_st
				do x = -this.dim_st, this.dim_st
					this.f_cor(x, y, face_idx)= 2*omega_cor*dsin(grid_points(2, x, y, face_idx))
				end do
			end do
		end do


		face_idx = 1
		do x = -dim, dim
			do y = -dim, dim

if(y+1 > dim)then ! y
call distance_sphere(r_sphere, grid_points(:, y, x, face_idx), grid_points(:, y-1, x, face_idx), dist)
this.h_dist(1, y, x) = dist
else
call distance_sphere(r_sphere, grid_points(:, y+1, x, face_idx), grid_points(:, y, x, face_idx), dist)
this.h_dist(1, y, x) = dist
end if

if(x+1 > dim)then ! x
call distance_sphere(r_sphere, grid_points(:, y, x, face_idx), grid_points(:, y, x-1, face_idx), dist)
this.h_dist(2, y, x) = dist
else
call distance_sphere(r_sphere, grid_points(:, y, x+1, face_idx), grid_points(:, y, x, face_idx), dist)
this.h_dist(2, y, x) = dist
end if

if(y-1 < -dim)then ! y
call distance_sphere(r_sphere, grid_points(:, y+1, x, face_idx), grid_points(:, y, x, face_idx), dist)
this.h_dist(3, y, x) = dist
else
call distance_sphere(r_sphere, grid_points(:, y, x, face_idx), grid_points(:, y-1, x, face_idx), dist)
this.h_dist(3, y, x) = dist
end if

if(x-1 < -dim)then ! x
call distance_sphere(r_sphere, grid_points(:, y, x+1, face_idx), grid_points(:, y, x, face_idx), dist)
this.h_dist(4, y, x) = dist
else
call distance_sphere(r_sphere, grid_points(:, y, x, face_idx), grid_points(:, y, x-1, face_idx), dist)
this.h_dist(4, y, x) = dist
end if

			end do
		end do

	end subroutine


	subroutine deinit(this)
		Class(g_var) :: this

		if (Allocated(this.h_dist)) Deallocate(this.h_dist)
		if (Allocated(this.f_cor)) Deallocate(this.f_cor)

	end subroutine



	real(8) function partial_c1_x(this, fun, x, y)
		Class(g_var) :: this
		real(8), intent(in) :: fun(1:3)
		integer, intent(in) :: x, y
		partial_c1_x = (fun(3) - fun(1))/(this.h_dist(4, y, x) + this.h_dist(2, y, x))

	end function



	real(8) function partial_c1_y(this, fun, x, y)
		Class(g_var) :: this
		real(8), intent(in) :: fun(1:3)
		integer, intent(in) :: x, y
		partial_c1_y = (fun(3) - fun(1))/(this.h_dist(3, y, x) + this.h_dist(1, y, x))

	end function


end module