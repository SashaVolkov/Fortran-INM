module special_variables

! use geometry

implicit none

	Private
	Public :: variables

	Type variables

		Real(8), Allocatable :: h_height(:, :, :)
		Real(8), Allocatable :: u_vel(:, :, :)
		Real(8), Allocatable :: v_vel(:, :, :)
		Real(8), Allocatable :: distance_grid(:, :, :, :)
		Real(8), Allocatable :: f_cor(:, :)
		Real(8), Allocatable :: alpha(:, :)
		Real(8), Allocatable :: beta(:, :)
		integer(4) dim, step

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Private :: alloc => alloc
		Procedure, Public :: deinit => deinit
		Procedure, Public :: eq => equal
		Procedure, Public :: start_conditions => start_conditions
	End Type


CONTAINS


	subroutine alloc(this)

		Class(variables) :: this

		Allocate(this.h_height(1:6, -this.dim:this.dim, -this.dim:this.dim))
		Allocate(this.u_vel(1:6, -this.dim:this.dim, -this.dim:this.dim))
		Allocate(this.v_vel(1:6, -this.dim:this.dim, -this.dim:this.dim))

		Allocate(this.distance_grid(1:6, -this.dim:this.dim, -this.dim:this.dim, 4*this.step))

		Allocate(this.f_cor(1:6, -this.dim:this.dim)) ! Only longitude
		Allocate(this.alpha(1:6, -this.dim:this.dim))
		Allocate(this.beta(1:6, -this.dim:this.dim))

	end subroutine


	subroutine init(this, grid_points, dim, step, omega_cor, r_sphere, g)

		Class(variables) :: this
		integer(4), intent(in) :: dim, step ! dimension
		real(8), intent(in) :: grid_points(1:6, -dim:dim, -dim:dim, 1:2), omega_cor, r_sphere, g

		integer(4) face_idx, x, y

		this.dim = dim+step
		this.step = step

		call this.alloc()

		do face_idx = 1, 6
			do x = -this.dim, this.dim ! Only longitude
				this.f_cor(face_idx, x)= 2*omega_cor*dsin(grid_points(face_idx, 0, x, 2))
				this.alpha(face_idx, x) = 1d0/(r_sphere*dcos(grid_points(face_idx, 0, x, 2)))
				this.beta(face_idx, x) = dtan(grid_points(face_idx, 0, x, 2))/r_sphere
			end do
		end do



		do face_idx = 1, 6
			do y = -dim, dim
				do x = -dim, dim

if(y+1 > dim)then ! latitude
this.distance_grid(face_idx, y, x, 1) = abs(grid_points(face_idx, y, x, 1)*r_sphere - grid_points(face_idx, y-1, x, 1)*r_sphere)
else
this.distance_grid(face_idx, y, x, 1) = abs(grid_points(face_idx, y+1, x, 1)*r_sphere - grid_points(face_idx, y, x, 1)*r_sphere)
end if

if(x+1 > dim)then ! longitude
this.distance_grid(face_idx, y, x, 2) = abs(grid_points(face_idx, y, x, 2)*r_sphere - grid_points(face_idx, y, x-1, 2)*r_sphere)
else
this.distance_grid(face_idx, y, x, 2) = abs(grid_points(face_idx, y, x+1, 2)*r_sphere - grid_points(face_idx, y, x, 2)*r_sphere)
end if

if(y-1 < -dim)then ! latitude
this.distance_grid(face_idx, y, x, 3) = abs(grid_points(face_idx, y+1, x, 1)*r_sphere - grid_points(face_idx, y, x, 1)*r_sphere)
else
this.distance_grid(face_idx, y, x, 3) = abs(grid_points(face_idx, y, x, 1)*r_sphere - grid_points(face_idx, y-1, x, 1)*r_sphere)
end if

if(x-1 < -dim)then ! longitude
this.distance_grid(face_idx, y, x, 4) = abs(grid_points(face_idx, y, x+1, 2)*r_sphere - grid_points(face_idx, y, x, 2)*r_sphere)
else
this.distance_grid(face_idx, y, x, 4) = abs(grid_points(face_idx, y, x, 2)*r_sphere - grid_points(face_idx, y, x-1, 2)*r_sphere)
end if

				end do
			end do
		end do

	print '("  Grid step real = ", f10.2, " m")', this.distance_grid(2, 300, 300, 4)


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
		integer(4), intent(in) :: dim ! dimension
		real(8), intent(in) :: grid_points(1:6, -dim:dim, -dim:dim, 1:2)

		integer(4) face_idx, x, y


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
				this.h_height(face_idx, y, x) = 0
			end do
		end do

	end subroutine



end module