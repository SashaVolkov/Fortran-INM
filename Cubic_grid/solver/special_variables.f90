module special_variables

use geometry

implicit none

	Private
	Public :: variables

	Type variables

		Real(8), Allocatable :: h_height(:, :, :)
		Real(8), Allocatable :: u_vel(:, :, :)
		Real(8), Allocatable :: v_vel(:, :, :)
		Real(8), Allocatable :: distance_grid(:, :, :, :)
		Real(8), Allocatable :: f_cor(:, :, :)
		Real(8), Allocatable :: alpha(:, :, :)
		Real(8), Allocatable :: beta(:, :, :)
		real(8) g, height
		integer(4) dim, step, dim_st

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Private :: alloc => alloc
		Procedure, Private :: const_def => const_def
		Procedure, Public :: deinit => deinit
	End Type


CONTAINS



	subroutine init(this, grid_points, dim, step, omega_cor, r_sphere, g, height)

		Class(variables) :: this
		integer(4), intent(in) :: dim, step ! dimension
		real(8), intent(in) :: grid_points(1:2, -dim:dim, -dim:dim, 1:6), omega_cor, r_sphere, g, height

		this.dim = dim;  this.step = step;  this.g = g;
		this.height = height;  this.dim_st = dim + step

		call this.alloc()
		call this.const_def(grid_points, dim, step, omega_cor, r_sphere)

	! print '("  Grid step real = ", f10.3, " m")', this.distance_grid(3, 300, 200, 2)

	end subroutine



		subroutine alloc(this)

			Class(variables) :: this

			Allocate(this.h_height(-this.dim_st:this.dim_st, -this.dim_st:this.dim_st, 1:6))
			Allocate(this.u_vel(-this.dim_st:this.dim_st, -this.dim_st:this.dim_st, 1:6))
			Allocate(this.v_vel(-this.dim_st:this.dim_st, -this.dim_st:this.dim_st, 1:6))

			Allocate(this.distance_grid(4*this.step, -this.dim_st:this.dim_st, -this.dim_st:this.dim_st, 1:6))

			Allocate(this.f_cor(-this.dim_st:this.dim_st, -this.dim_st:this.dim_st, 1:6))
			Allocate(this.alpha(-this.dim_st:this.dim_st, -this.dim_st:this.dim_st, 1:6))
			Allocate(this.beta(-this.dim_st:this.dim_st, -this.dim_st:this.dim_st, 1:6))

		end subroutine



		subroutine const_def(this, grid_points, dim, step, omega_cor, r_sphere)
			Class(variables) :: this
			integer(4), intent(in) :: dim, step ! dimension
			real(8), intent(in) :: grid_points(1:2, -dim:dim, -dim:dim, 1:6), omega_cor, r_sphere
			real(8) dist
			integer(4) face_idx, x, y


			do face_idx = 1, 6 ! Only longitude
				do y = -this.dim_st, this.dim_st
					do x = -this.dim_st, this.dim_st
						this.f_cor(x, y, face_idx)= 2*omega_cor*dsin(grid_points(2, x, y, face_idx))
						this.alpha(x, y, face_idx) = 1d0/(r_sphere*dcos(grid_points(2, x, y, face_idx)))
						this.beta(x, y, face_idx) = dtan(grid_points(2, x, y, face_idx))/r_sphere
					end do
				end do
			end do


		do face_idx = 1, 6 ! calculate all distances between points
			do x = -dim, dim
				do y = -dim, dim

if(y+1 > dim)then ! latitude
call distance_sphere(r_sphere, grid_points(:, y, x, face_idx), grid_points(:, y-1, x, face_idx), dist)
this.distance_grid(1, y, x, face_idx) = dist
else
call distance_sphere(r_sphere, grid_points(:, y+1, x, face_idx), grid_points(:, y, x, face_idx), dist)
this.distance_grid(1, y, x, face_idx) = dist
end if

if(x+1 > dim)then ! longitude
call distance_sphere(r_sphere, grid_points(:, y, x, face_idx), grid_points(:, y, x-1, face_idx), dist)
this.distance_grid(2, y, x, face_idx) = dist
else
call distance_sphere(r_sphere, grid_points(:, y, x+1, face_idx), grid_points(:, y, x, face_idx), dist)
this.distance_grid(2, y, x, face_idx) = dist
end if

if(y-1 < -dim)then ! latitude
call distance_sphere(r_sphere, grid_points(:, y+1, x, face_idx), grid_points(:, y, x, face_idx), dist)
this.distance_grid(3, y, x, face_idx) = dist
else
call distance_sphere(r_sphere, grid_points(:, y, x, face_idx), grid_points(:, y-1, x, face_idx), dist)
this.distance_grid(3, y, x, face_idx) = dist
end if

if(x-1 < -dim)then ! longitude
call distance_sphere(r_sphere, grid_points(:, y, x+1, face_idx), grid_points(:, y, x, face_idx), dist)
this.distance_grid(4, y, x, face_idx) = dist
else
call distance_sphere(r_sphere, grid_points(:, y, x, face_idx), grid_points(:, y, x-1, face_idx), dist)
this.distance_grid(4, y, x, face_idx) = dist
end if

				end do
			end do
		end do

		end subroutine


	subroutine deinit(this)
		Class(variables) :: this

		if (Allocated(this.h_height)) Deallocate(this.h_height)
		if (Allocated(this.u_vel)) Deallocate(this.u_vel)
		if (Allocated(this.v_vel)) Deallocate(this.v_vel)
		if (Allocated(this.f_cor)) Deallocate(this.f_cor)
		if (Allocated(this.alpha)) Deallocate(this.alpha)
		if (Allocated(this.beta)) Deallocate(this.beta)

	end subroutine

end module