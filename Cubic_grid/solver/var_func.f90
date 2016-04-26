module var_func

 use special_variables, Only: variables

implicit none

	Private
	Public :: func

	Type func
		CONTAINS
		Procedure, Public :: equal => equal
		Procedure, Public :: start_conditions => start_conditions
		Procedure, Public :: borders => borders
		Procedure, Public :: corner_zero => corner_zero
	End Type


CONTAINS



	subroutine equal(f, var_pr, var)

		Class(func) :: f
		Class(variables) :: var_pr, var
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



	subroutine start_conditions(f, var_pr)

		Class(func) :: f
		Class(variables) :: var_pr
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
				var_pr.h_height(x, y, 4) = h0*exp(-((((10.0/dim)*(x*0.5))**2)+(((10.0/dim)*(y*0.5))**2)))
				! var_pr.h_height(x, y, 2) = 0
			end do
		end do

	end subroutine



	subroutine borders(f, var, var_grey)
	Class(func) :: f
	Class(variables) :: var, var_grey
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



	subroutine corner_zero(f, var)
	Class(func) :: f
	Class(variables) :: var
	integer(4) dim, step, dim_st, face_idx

	dim = var.dim;  step = var.step;  dim_st = var.dim_st

	do face_idx = 1,6
		var.h_height(dim_st, dim_st, face_idx) = 0
		var.u_vel(dim_st, dim_st, face_idx) = 0
		var.v_vel(dim_st, dim_st, face_idx) = 0

		var.h_height(-dim_st, dim_st, face_idx) = 0
		var.u_vel(-dim_st, dim_st, face_idx) = 0
		var.v_vel(-dim_st, dim_st, face_idx) = 0

		var.h_height(dim_st, -dim_st, face_idx) = 0
		var.u_vel(dim_st, -dim_st, face_idx) = 0
		var.v_vel(dim_st, -dim_st, face_idx) = 0

		var.h_height(-dim_st, -dim_st, face_idx) = 0
		var.u_vel(-dim_st, -dim_st, face_idx) = 0
		var.v_vel(-dim_st, -dim_st, face_idx) = 0
	end do

	end subroutine



end module