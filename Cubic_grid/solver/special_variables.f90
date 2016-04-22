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
		integer(4) dim, step

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Private :: alloc => alloc
		Procedure, Private :: const_def => const_def
		Procedure, Public :: deinit => deinit
		Procedure, Public :: eq => equal
		Procedure, Public :: start_conditions => start_conditions
	End Type


CONTAINS



	subroutine init(this, grid_points, dim, step, omega_cor, r_sphere, g, height)

		Class(variables) :: this
		integer(4), intent(in) :: dim, step ! dimension
		real(8), intent(in) :: grid_points(1:6, -dim:dim, -dim:dim, 1:2), omega_cor, r_sphere, g, height

		this.dim = dim+step;  this.step = step;  this.g = g;  this.height = height

		call this.alloc()
		call this.const_def(grid_points, dim, step, omega_cor, r_sphere)

	print '("  Grid step real = ", f10.3, " m")', this.distance_grid(2, 200, 200, 3)

	end subroutine



		subroutine alloc(this)

			Class(variables) :: this

			Allocate(this.h_height(1:6, -this.dim:this.dim, -this.dim:this.dim))
			Allocate(this.u_vel(1:6, -this.dim:this.dim, -this.dim:this.dim))
			Allocate(this.v_vel(1:6, -this.dim:this.dim, -this.dim:this.dim))

			Allocate(this.distance_grid(1:6, -this.dim:this.dim, -this.dim:this.dim, 4*this.step))

			Allocate(this.f_cor(1:6, -this.dim:this.dim, -this.dim:this.dim))
			Allocate(this.alpha(1:6, -this.dim:this.dim, -this.dim:this.dim))
			Allocate(this.beta(1:6, -this.dim:this.dim, -this.dim:this.dim))

		end subroutine



		subroutine const_def(this, grid_points, dim, step, omega_cor, r_sphere)
			Class(variables) :: this
			integer(4), intent(in) :: dim, step ! dimension
			real(8), intent(in) :: grid_points(1:6, -dim:dim, -dim:dim, 1:2), omega_cor, r_sphere
			real(8) dist
			integer(4) face_idx, x, y


			do face_idx = 1, 6 ! Only longitude
				do y = -this.dim, this.dim
					do x = -this.dim, this.dim
						this.f_cor(face_idx, y, x)= 2*omega_cor*dsin(grid_points(face_idx, y, x, 2))
						this.alpha(face_idx, y, x) = 1d0/(r_sphere*dcos(grid_points(face_idx, y, x, 2)))
						this.beta(face_idx, y, x) = dtan(grid_points(face_idx, y, x, 2))/r_sphere
					end do
				end do
			end do


		do face_idx = 1, 6 ! calculate all distances between points
			do y = -dim, dim
				do x = -dim, dim

if(y+1 > dim)then ! latitude
call distance_sphere(r_sphere, grid_points(face_idx, y, x, :), grid_points(face_idx, y-1, x, :), dist)
this.distance_grid(face_idx, y, x, 1) = dist
else
call distance_sphere(r_sphere, grid_points(face_idx, y+1, x, :), grid_points(face_idx, y, x, :), dist)
this.distance_grid(face_idx, y, x, 1) = dist
end if

if(x+1 > dim)then ! longitude
call distance_sphere(r_sphere, grid_points(face_idx, y, x, :), grid_points(face_idx, y, x-1, :), dist)
this.distance_grid(face_idx, y, x, 2) = dist
else
call distance_sphere(r_sphere, grid_points(face_idx, y, x+1, :), grid_points(face_idx, y, x, :), dist)
this.distance_grid(face_idx, y, x, 2) = dist
end if

if(y-1 < -dim)then ! latitude
call distance_sphere(r_sphere, grid_points(face_idx, y+1, x, :), grid_points(face_idx, y, x, :), dist)
this.distance_grid(face_idx, y, x, 3) = dist
else
call distance_sphere(r_sphere, grid_points(face_idx, y, x, :), grid_points(face_idx, y-1, x, :), dist)
this.distance_grid(face_idx, y, x, 3) = dist
end if

if(x-1 < -dim)then ! longitude
call distance_sphere(r_sphere, grid_points(face_idx, y, x+1, :), grid_points(face_idx, y, x, :), dist)
this.distance_grid(face_idx, y, x, 4) = dist
else
call distance_sphere(r_sphere, grid_points(face_idx, y, x, :), grid_points(face_idx, y, x-1, :), dist)
this.distance_grid(face_idx, y, x, 4) = dist
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



	subroutine equal(this, that)

		Class(variables) :: this, that
		integer(4) dim, face_idx, x, y
		dim = this.dim

		do face_idx = 1, 6
			do y = -dim, dim
				do x = -dim, dim
				this.h_height(face_idx, y, x)=that.h_height(face_idx, y, x)
				this.u_vel(face_idx, y, x)=that.u_vel(face_idx, y, x)
				this.v_vel(face_idx, y, x)=that.v_vel(face_idx, y, x)
				end do
			end do
		end do

	end subroutine



	subroutine start_conditions(this, grid_points, dim)

		Class(variables) :: this
		integer(4), intent(in) :: dim
		real(8), intent(in) :: grid_points(1:6, -dim:dim, -dim:dim, 1:2)
		real(8) h0

		integer(4) face_idx, x, y

		h0 = this.height

		do face_idx = 1, 6
			do y = -dim, dim
				do x = -dim, dim
					this.h_height(face_idx, y, x) = 0
					this.u_vel(face_idx, y, x) = 0
					this.v_vel(face_idx, y, x) = 0
				end do
			end do
		end do

		face_idx = 2
		do y = -dim, dim
			do x = -dim, dim
this.h_height(face_idx, y, x) = h0*exp(-((((10.0/dim)*(x-dim*0.5))**2)+(((10.0/dim)*(y-dim*0.5))**2)))
			end do
		end do

	end subroutine



end module