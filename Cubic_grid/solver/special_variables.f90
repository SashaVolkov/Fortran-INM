module special_variables

! use 

implicit none

	Private
	Public :: variables

	Type variables

		Real(8), Allocatable :: h_height(:, :, :)
		Real(8), Allocatable :: u_vel(:, :, :)
		Real(8), Allocatable :: v_vel(:, :, :)
		Real(8), Allocatable :: f_cor(:, :, :)
		Real(8), Allocatable :: alpha(:, :, :)
		Real(8), Allocatable :: beta(:, :, :)
		integer(4) dim

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Public :: deinit => deinit
		Procedure, Public :: eq => equal
	End Type


CONTAINS

	subroutine init(this, grid_points, dim, omega_cor, r_sphere, g)

		Class(variables) :: this
		integer(4), intent(in) :: dim ! dimension
		real(8), intent(in) :: grid_points(1:6, -dim:dim, -dim:dim, 1:2), omega_cor, r_sphere, g

		integer(4) face_idx, x, y

		this.dim = dim

		Allocate(this.h_height(1:6, -dim:dim, -dim:dim))
		Allocate(this.u_vel(1:6, -dim:dim, -dim:dim))
		Allocate(this.v_vel(1:6, -dim:dim, -dim:dim))
		Allocate(this.f_cor(1:6, -dim:dim, -dim:dim))
		Allocate(this.alpha(1:6, -dim:dim, -dim:dim))
		Allocate(this.beta(1:6, -dim:dim, -dim:dim))


		do face_idx = 1, 6
			do y = -dim, dim
				do x = -dim, dim
				this.f_cor(face_idx, y, x)= 2*omega_cor*dsin(grid_points(face_idx, y, x, 1))
				this.alpha(face_idx, y, x) = 1d0/(r_sphere*dcos(grid_points(face_idx, y, x, 1)))
				this.beta(face_idx, y, x) = dtan(grid_points(face_idx, y, x, 1))/r_sphere
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



	subroutine start_conditions(this, grid_points, dim, omega_cor)

		Class(variables) :: this
		integer(4), intent(in) :: dim ! dimension
		real(8), intent(in) :: grid_points(1:6, -dim:dim, -dim:dim, 1:2), omega_cor

		integer(4) face_idx, x, y


		do face_idx = 1, 6
			do y = -dim, dim
				do x = -dim, dim



				end do
			end do
		end do

	end subroutine



end module