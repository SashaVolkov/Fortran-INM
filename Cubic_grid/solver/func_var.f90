module func_var

use geometry

implicit none

	Private
	Public :: f_var

	Type f_var

		Real(8), Allocatable :: h_height(:, :, :)
		Real(8), Allocatable :: u_vel(:, :, :)
		Real(8), Allocatable :: v_vel(:, :, :)
		! Real(8), Allocatable :: distance_grid(:, :, :, :)
		real(8) g, height
		integer(4) dim, step, dim_st

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Private :: alloc => alloc
		Procedure, Public :: deinit => deinit
		Procedure, Public :: equal => equal
		Procedure, Public :: start_conditions => start_conditions
		Procedure, Public :: borders => borders
		Procedure, Public :: corner_zero => corner_zero
	End Type


CONTAINS



	subroutine init(this, grid_points, dim, step, omega_cor, r_sphere, g, height)

		Class(f_var) :: this
		integer(4), intent(in) :: dim, step ! dimension
		real(8), intent(in) :: grid_points(1:2, -dim:dim, -dim:dim, 1:6), omega_cor, r_sphere, g, height

		this.dim = dim;  this.step = step;  this.g = g;
		this.height = height;  this.dim_st = dim + step

		call this.alloc()
		! call this.const_def(grid_points, dim, step, omega_cor, r_sphere)

	! print '("  Grid step real = ", f10.3, " m")', this.distance_grid(3, 300, 200, 2)

	end subroutine



		subroutine alloc(this)

			Class(f_var) :: this

			Allocate(this.h_height(-this.dim_st:this.dim_st, -this.dim_st:this.dim_st, 1:6))
			Allocate(this.u_vel(-this.dim_st:this.dim_st, -this.dim_st:this.dim_st, 1:6))
			Allocate(this.v_vel(-this.dim_st:this.dim_st, -this.dim_st:this.dim_st, 1:6))

			! Allocate(this.distance_grid(4*this.step, -this.dim_st:this.dim_st, -this.dim_st:this.dim_st, 1:6))

		end subroutine



	subroutine deinit(this)
		Class(f_var) :: this

		! if (Allocated(this.distance_grid)) Deallocate(this.distance_grid)
		if (Allocated(this.h_height)) Deallocate(this.h_height)
		if (Allocated(this.u_vel)) Deallocate(this.u_vel)
		if (Allocated(this.v_vel)) Deallocate(this.v_vel)

	end subroutine



		subroutine equal(var_pr, var)

			Class(f_var) :: var_pr, var
			integer(4) dim, face_idx, x, y
			dim = var_pr.dim_st

			do face_idx = 1, 6
				do y = -dim, dim
					do x = -dim, dim
					var_pr.h_height(x, y, face_idx)=var.h_height(x, y, face_idx)
					var_pr.u_vel(x, y, face_idx)=var.u_vel(x, y, face_idx)
					var_pr.v_vel(x, y, face_idx)=var.v_vel(x, y, face_idx)
					end do
				end do
			end do

		end subroutine



		subroutine start_conditions(var_pr)

			Class(f_var) :: var_pr
			integer(4) dim
			real(8) h0

			integer(4) face_idx, x, y

			h0 = var_pr.height;  dim = var_pr.dim_st

			do face_idx = 1, 6
				do y = -dim, dim
					do x = -dim, dim
						var_pr.h_height(x, y, face_idx) = 0
						var_pr.u_vel(x, y, face_idx) = 0
						var_pr.v_vel(x, y, face_idx) = 0
					end do
				end do
			end do


			do y = -dim, dim
				do x = -dim, dim
					var_pr.h_height(x, y, 2) = h0*exp(-((((10.0/dim)*(x*0.5))**2)+(((10.0/dim)*(y*0.5))**2)))
	! 				var_pr.h_height(x, y, 4) = h0*exp(-((((10.0/dim)*(x*0.5))**2)+(((10.0/dim)*(y*0.5))**2)))
					! var_pr.h_height(x, y, 2) = 0
				end do
			end do

		end subroutine



		subroutine borders(var, var_grey)

		Class(f_var) :: var, var_grey
		integer(4) dim, step, dim_st, x, y, i, j

		dim = var.dim;  step = var.step;  dim_st = var.dim_st


		do y = 0, step
			do x = -dim, dim

				var_grey.h_height(x, -dim - y, 1) = var.h_height(x, dim - y, 2)
				var_grey.u_vel(x, -dim - y, 1) = var.u_vel(x, dim - y, 2)
				var_grey.v_vel(x, -dim - y, 1) = var.v_vel(x, dim - y, 2)

				var_grey.h_height(x, dim + y, 2) = var.h_height(x, -dim + y, 1)
				var_grey.u_vel(x, dim + y, 2) = var.u_vel(x, -dim + y, 1)
				var_grey.v_vel(x, dim + y, 2) = var.v_vel(x, -dim + y, 1)


				var_grey.h_height(x, dim + y, 6) = var.h_height(x, -dim + y, 2)
				var_grey.u_vel(x, dim + y, 6) = var.u_vel(x, -dim + y, 2)
				var_grey.v_vel(x, dim + y, 6) = var.v_vel(x, -dim + y, 2)

				var_grey.h_height(x, -dim - y, 2) = var.h_height(x, dim - y, 6)
				var_grey.u_vel(x, -dim - y, 2) = var.u_vel(x, dim - y, 6)
				var_grey.v_vel(x, -dim - y, 2) = var.v_vel(x, dim - y, 6)

			end do
		end do


		do y = -dim, dim
			do x = 0, step

				var_grey.h_height(dim - x, y, 5) = var.h_height(-dim + x, y, 2)
				var_grey.u_vel(dim - x, y, 5) = var.u_vel(-dim + x, y, 2)
				var_grey.v_vel(dim - x, y, 5) = var.v_vel(-dim + x, y, 2)

				var_grey.h_height(-dim - x, y, 2) = var.h_height(dim + x, y, 5)
				var_grey.u_vel(-dim - x, y, 2) = var.u_vel(dim + x, y, 5)
				var_grey.v_vel(-dim - x, y, 2) = var.v_vel(dim + x, y, 5)


				var_grey.h_height(-dim - x, y, 3) = var.h_height(dim - x, y, 2)
				var_grey.u_vel(-dim - x, y, 3) = var.u_vel(dim - x, y, 2)
				var_grey.v_vel(-dim - x, y, 3) = var.v_vel(dim - x, y, 2)

				var_grey.h_height(dim + x, y, 2) = var.h_height(-dim + x, y, 3)
				var_grey.u_vel(dim + x, y, 2) = var.u_vel(-dim + x, y, 3)
				var_grey.v_vel(dim + x, y, 2) = var.v_vel(-dim + x, y, 3)

			end do
		end do


		do y = 0, step
			do x = -dim, dim

			var_grey.h_height(-x, dim + y, 1) = var.h_height(x, dim - y, 4)
			var_grey.u_vel(-x, dim + y, 1) = -var.u_vel(x, dim - y, 4)
			var_grey.v_vel(-x, dim + y, 1) = -var.v_vel(x, dim - y, 4)

				var_grey.h_height(x, dim + y, 4) = var.h_height(-x, dim - y, 1)
				var_grey.u_vel(x, dim + y, 4) = -var.u_vel(-x, dim - y, 1)
				var_grey.v_vel(x, dim + y, 4) = -var.v_vel(-x, dim - y, 1)


				var_grey.h_height(-x, -dim - y, 6) = var.h_height(x, -dim + y, 4)
				var_grey.u_vel(-x, -dim - y, 6) = -var.u_vel(x, -dim + y, 4)
				var_grey.v_vel(-x, -dim - y, 6) = -var.v_vel(x, -dim + y, 4)

				var_grey.h_height(x, -dim - y, 4) = var.h_height(-x, -dim + y, 6)
				var_grey.u_vel(x, -dim - y, 4) = -var.u_vel(-x, -dim + y, 6)
				var_grey.v_vel(x, -dim - y, 4) = -var.v_vel(-x, -dim + y, 6)

			end do
		end do


		do y = -dim, dim
			do x = 0, step

				var_grey.h_height(dim + x, y, 3) = var.h_height(-dim + x, y, 4)
				var_grey.u_vel(dim + x, y, 3) = var.u_vel(-dim + x, y, 4)
				var_grey.v_vel(dim + x, y, 3) = var.v_vel(-dim + x, y, 4)

				var_grey.h_height(-dim - x, y, 4) = var.h_height(dim - x, y, 3)
				var_grey.u_vel(-dim - x, y, 4) = var.u_vel(dim - x, y, 3)
				var_grey.v_vel(-dim - x, y, 4) = var.v_vel(dim - x, y, 3)


				var_grey.h_height(-dim - x, y, 5) = var.h_height(dim - x, y, 4)
				var_grey.u_vel(-dim - x, y, 5) = var.u_vel(dim - x, y, 4)
				var_grey.v_vel(-dim - x, y, 5) = var.v_vel(dim - x, y, 4)

				var_grey.h_height(dim + x, y, 4) = var.h_height(-dim + x, y, 5)
				var_grey.u_vel(dim + x, y, 4) = var.u_vel(-dim + x, y, 5)
				var_grey.v_vel(dim + x, y, 4) = var.v_vel(-dim + x, y, 5)

			end do
		end do


		do i = -dim, dim
			do j = 0, step

				var_grey.h_height(-dim - j, -i, 1) = var.h_height(i, dim - j, 5)
				var_grey.u_vel(-dim - j, -i, 1) = var.v_vel(i, dim - j, 5)
				var_grey.v_vel(-dim - j, -i, 1) = -var.u_vel(i, dim - j, 5)

				var_grey.h_height(i, dim + j, 5) = var.h_height(-dim + j, -i, 1)
				var_grey.u_vel(i, dim + j, 5) = -var.v_vel(-dim + j, -i, 1)
				var_grey.v_vel(i, dim + j, 5) = var.u_vel(-dim + j, -i, 1)


				var_grey.h_height(-dim - j, i, 6) = var.h_height(i, -dim + j, 5)
				var_grey.u_vel(-dim - j, i, 6) = -var.v_vel(i, -dim + j, 5)
				var_grey.v_vel(-dim - j, i, 6) = var.u_vel(i, -dim + j, 5)

				var_grey.h_height(i, -dim - j, 5) = var.h_height(-dim + j, i, 6)
				var_grey.u_vel(i, -dim - j, 5) = var.v_vel(-dim + j, i, 6)
				var_grey.v_vel(i, -dim - j, 5) = -var.u_vel(-dim + j, i, 6)



				var_grey.h_height(dim + j, i, 1) = var.h_height(i, dim - j, 3)
				var_grey.u_vel(dim + j, i, 1) = -var.v_vel(i, dim - j, 3)
				var_grey.v_vel(dim + j, i, 1) = var.u_vel(i, dim - j, 3)

				var_grey.h_height(i, dim + j, 3) = var.h_height(dim - j, i, 1)
				var_grey.u_vel(i, dim + j, 3) = var.v_vel(dim - j, i, 1)
				var_grey.v_vel(i, dim + j, 3) = -var.u_vel(dim - j, i, 1)


				var_grey.h_height(dim + j, -i, 6) = var.h_height(i, -dim + j, 3)
				var_grey.u_vel(dim + j, -i, 6) = var.v_vel(i, -dim + j, 3)
				var_grey.v_vel(dim + j, -i, 6) = -var.u_vel(i, -dim + j, 3)

				var_grey.h_height(i, -dim - j, 3) = var.h_height(dim - j, -i, 6)
				var_grey.u_vel(i, -dim - j, 3) = -var.v_vel(dim - j, -i, 6)
				var_grey.v_vel(i, -dim - j, 3) = var.u_vel(dim - j, -i, 6)

			end do
		end do


		end subroutine



		subroutine corner_zero(var)

		Class(f_var) :: var
		integer(4) dim_st, face_idx, x,y

		dim_st = var.dim_st

		do face_idx = 1,6
			do y = -dim_st, dim_st, 2*dim_st+1
				do x = -dim_st, dim_st, 2*dim_st+1
			var.h_height(x, y, face_idx) = 0
			var.u_vel(x, y, face_idx) = 0
			var.v_vel(x, y, face_idx) = 0

				end do
			end do
		end do

		end subroutine



end module