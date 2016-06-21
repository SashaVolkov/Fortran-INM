module func_var

	use parallel_cubic, Only: parallel

implicit none

	Private
	Public :: f_var

	Type f_var

		Real(8), Allocatable :: h_height(:, :, :)
		Real(8), Allocatable :: x_vel(:, :, :)
		Real(8), Allocatable :: y_vel(:, :, :)
		! Real(8), Allocatable :: distance_grid(:, :, :, :)
		real(8) height
		integer(4) step, dim, Xsize, Ysize
		integer(4) ns_x, ns_y, nf_x, nf_y, first_x, first_y, last_x, last_y

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Private :: alloc => alloc
		Procedure, Public :: deinit => deinit
		Procedure, Public :: equal => equal
		Procedure, Public :: start_conditions => start_conditions
		! Procedure, Public :: borders => borders
	End Type


CONTAINS



	subroutine init(this, paral, step, height)

		Class(f_var) :: this
		Class(parallel) :: paral
		integer(4), intent(in) :: step
		real(8), intent(in) :: height


		this.ns_x = paral.ns_xy(1);  this.ns_y = paral.ns_xy(2)
		this.nf_x = paral.nf_xy(1);  this.nf_y = paral.nf_xy(2)

		this.first_x = paral.ns_xy(1) - step;  this.first_y = paral.ns_xy(2) - step
		this.last_x = paral.nf_xy(1) + step;  this.last_y = paral.nf_xy(2) + step

		this.Xsize = paral.Xsize;  this.Ysize = paral.Ysize

		this.step = step;  this.height = height;  this.dim = paral.dim


		call this.alloc()

	end subroutine



	subroutine alloc(this)

		Class(f_var) :: this

		Allocate(this.h_height(this.first_x:this.last_x, this.first_y:this.last_y, 6))
		Allocate(this.x_vel(this.first_x:this.last_x, this.first_y:this.last_y, 6))
		Allocate(this.y_vel(this.first_x:this.last_x, this.first_y:this.last_y, 6))

	end subroutine



	subroutine deinit(this)
		Class(f_var) :: this

		if (Allocated(this.h_height)) Deallocate(this.h_height)
		if (Allocated(this.x_vel)) Deallocate(this.x_vel)
		if (Allocated(this.y_vel)) Deallocate(this.y_vel)

	end subroutine



	subroutine equal(var_pr, var)

		Class(f_var) :: var_pr, var
		integer(4) x, y

			! do y = var_pr.first_y, var_pr.last_y
			! 	do x = var_pr.first_x, var_pr.last_x
				var_pr.h_height(:, :, :)=var.h_height(:, :, :)
				var_pr.x_vel(:, :, :)=var.x_vel(:, :, :)
				var_pr.y_vel(:, :, :)=var.y_vel(:, :, :)
			! 	end do
			! end do


	end subroutine



	subroutine start_conditions(this)

		Class(f_var) :: this
		integer(4) dim
		real(8) h0

		integer(4) x, y, face

		h0 = this.height;  dim = this.dim

		do face = 1, 6

			do y = this.first_y, this.last_y
				do x = this.first_x, this.last_x
					this.h_height(x, y, face) = 0
					this.x_vel(x, y, face) = 0
					this.y_vel(x, y, face) = 0
				end do
			end do

			if ( face == 2 ) then
			do y = this.first_y, this.last_y
				do x = this.first_x, this.last_x
					this.h_height(x, y, face) =&
					 h0*exp(-((((10.0/dim)*((x-dim)*0.5))**2)+(((10.0/dim)*((y-dim)*0.5))**2)))
				end do
			end do
			end if

			! if ( face == 1 ) then
			! do y = this.first_y, this.last_y
			! 	do x = this.first_x, this.last_x
			! 		this.h_height(x, y, face) = 1.0
			! 	end do
			! end do
			! end if



		end do


	end subroutine



		! integer(4) function borders(this, var, var_grey)

		! Class(f_var) :: this
		! Class(f_var) :: var(1:6), var_grey(1:6)
		! integer(4) dim, step, dim_st, x, y, i, j

		! dim = var(1).dim;  step = var(1).step;  dim_st = var(1).dim_st


		! do y = 0, step
		! 	do x = -dim, dim

		! 		var_grey(1).h_height(x, -dim - y) = var(2).h_height(x, dim - y)
		! 		var_grey(1).x_vel(x, -dim - y) = var(2).x_vel(x, dim - y)
		! 		var_grey(1).y_vel(x, -dim - y) = var(2).y_vel(x, dim - y)

		! 		var_grey(2).h_height(x, dim + y) = var(1).h_height(x, -dim + y)
		! 		var_grey(2).x_vel(x, dim + y) = var(1).x_vel(x, -dim + y)
		! 		var_grey(2).y_vel(x, dim + y) = var(1).y_vel(x, -dim + y)


		! 		var_grey(6).h_height(x, dim + y) = var(2).h_height(x, -dim + y)
		! 		var_grey(6).x_vel(x, dim + y) = var(2).x_vel(x, -dim + y)
		! 		var_grey(6).y_vel(x, dim + y) = var(2).y_vel(x, -dim + y)

		! 		var_grey(2).h_height(x, -dim - y) = var(6).h_height(x, dim - y)
		! 		var_grey(2).x_vel(x, -dim - y) = var(6).x_vel(x, dim - y)
		! 		var_grey(2).y_vel(x, -dim - y) = var(6).y_vel(x, dim - y)

		! 	end do
		! end do


		! do y = -dim, dim
		! 	do x = 0, step

		! 		var_grey(5).h_height(dim + x, y) = var(2).h_height(-dim + x, y)
		! 		var_grey(5).x_vel(dim + x, y) = var(2).x_vel(-dim + x, y)
		! 		var_grey(5).y_vel(dim + x, y) = var(2).y_vel(-dim + x, y)

		! 		var_grey(2).h_height(-dim - x, y) = var(5).h_height(dim - x, y)
		! 		var_grey(2).x_vel(-dim - x, y) = var(5).x_vel(dim - x, y)
		! 		var_grey(2).y_vel(-dim - x, y) = var(5).y_vel(dim - x, y)


		! 		var_grey(3).h_height(-dim - x, y) = var(2).h_height(dim - x, y)
		! 		var_grey(3).x_vel(-dim - x, y) = var(2).x_vel(dim - x, y)
		! 		var_grey(3).y_vel(-dim - x, y) = var(2).y_vel(dim - x, y)

		! 		var_grey(2).h_height(dim + x, y) = var(3).h_height(-dim + x, y)
		! 		var_grey(2).x_vel(dim + x, y) = var(3).x_vel(-dim + x, y)
		! 		var_grey(2).y_vel(dim + x, y) = var(3).y_vel(-dim + x, y)

		! 	end do
		! end do


		! do y = 0, step
		! 	do x = -dim, dim

		! 	var_grey(1).h_height(-x, dim + y) = var(4).h_height(x, dim - y)
		! 	var_grey(1).x_vel(-x, dim + y) = -var(4).x_vel(x, dim - y)
		! 	var_grey(1).y_vel(-x, dim + y) = -var(4).y_vel(x, dim - y)

		! 		var_grey(4).h_height(x, dim + y) = var(1).h_height(-x, dim - y)
		! 		var_grey(4).x_vel(x, dim + y) = -var(1).x_vel(-x, dim - y)
		! 		var_grey(4).y_vel(x, dim + y) = -var(1).y_vel(-x, dim - y)


		! 		var_grey(6).h_height(-x, -dim - y) = var(4).h_height(x, -dim + y)
		! 		var_grey(6).x_vel(-x, -dim - y) = -var(4).x_vel(x, -dim + y)
		! 		var_grey(6).y_vel(-x, -dim - y) = -var(4).y_vel(x, -dim + y)

		! 		var_grey(4).h_height(x, -dim - y) = var(6).h_height(-x, -dim + y)
		! 		var_grey(4).x_vel(x, -dim - y) = -var(6).x_vel(-x, -dim + y)
		! 		var_grey(4).y_vel(x, -dim - y) = -var(6).y_vel(-x, -dim + y)

		! 	end do
		! end do


		! do y = -dim, dim
		! 	do x = 0, step

		! 		var_grey(3).h_height(dim + x, y) = var(4).h_height(-dim + x, y)
		! 		var_grey(3).x_vel(dim + x, y) = var(4).x_vel(-dim + x, y)
		! 		var_grey(3).y_vel(dim + x, y) = var(4).y_vel(-dim + x, y)

		! 		var_grey(4).h_height(-dim - x, y) = var(3).h_height(dim - x, y)
		! 		var_grey(4).x_vel(-dim - x, y) = var(3).x_vel(dim - x, y)
		! 		var_grey(4).y_vel(-dim - x, y) = var(3).y_vel(dim - x, y)


		! 		var_grey(4).h_height(dim + x, y) = var(5).h_height(-dim + x, y)
		! 		var_grey(4).x_vel(dim + x, y) = var(5).x_vel(-dim + x, y)
		! 		var_grey(4).y_vel(dim + x, y) = var(5).y_vel(-dim + x, y)

		! 		var_grey(5).h_height(-dim - x, y) = var(4).h_height(dim - x, y)
		! 		var_grey(5).x_vel(-dim - x, y) = var(4).x_vel(dim - x, y)
		! 		var_grey(5).y_vel(-dim - x, y) = var(4).y_vel(dim - x, y)

		! 	end do
		! end do


		! do i = -dim, dim
		! 	do j = 0, step

		! 		var_grey(1).h_height(-dim - j, -i) = var(5).h_height(i, dim - j)
		! 		var_grey(1).x_vel(-dim - j, -i) = var(5).y_vel(i, dim - j)
		! 		var_grey(1).y_vel(-dim - j, -i) = -var(5).x_vel(i, dim - j)

		! 		var_grey(5).h_height(i, dim + j) = var(1).h_height(-dim + j, -i)
		! 		var_grey(5).x_vel(i, dim + j) = -var(1).y_vel(-dim + j, -i)
		! 		var_grey(5).y_vel(i, dim + j) = var(1).x_vel(-dim + j, -i)


		! 		var_grey(6).h_height(-dim - j, i) = var(5).h_height(i, -dim + j)
		! 		var_grey(6).x_vel(-dim - j, i) = -var(5).y_vel(i, -dim + j)
		! 		var_grey(6).y_vel(-dim - j, i) = var(5).x_vel(i, -dim + j)

		! 		var_grey(5).h_height(i, -dim - j) = var(6).h_height(-dim + j, i)
		! 		var_grey(5).x_vel(i, -dim - j) = var(6).y_vel(-dim + j, i)
		! 		var_grey(5).y_vel(i, -dim - j) = -var(6).x_vel(-dim + j, i)



		! 		var_grey(1).h_height(dim + j, i) = var(3).h_height(i, dim - j)
		! 		var_grey(1).x_vel(dim + j, i) = -var(3).y_vel(i, dim - j)
		! 		var_grey(1).y_vel(dim + j, i) = var(3).x_vel(i, dim - j)

		! 		var_grey(3).h_height(i, dim + j) = var(1).h_height(dim - j, i)
		! 		var_grey(3).x_vel(i, dim + j) = var(1).y_vel(dim - j, i)
		! 		var_grey(3).y_vel(i, dim + j) = -var(1).x_vel(dim - j, i)


		! 		var_grey(6).h_height(dim + j, -i) = var(3).h_height(i, -dim + j)
		! 		var_grey(6).x_vel(dim + j, -i) = var(3).y_vel(i, -dim + j)
		! 		var_grey(6).y_vel(dim + j, -i) = -var(3).x_vel(i, -dim + j)

		! 		var_grey(3).h_height(i, -dim - j) = var(6).h_height(dim - j, -i)
		! 		var_grey(3).x_vel(i, -dim - j) = -var(6).y_vel(dim - j, -i)
		! 		var_grey(3).y_vel(i, -dim - j) = var(6).x_vel(dim - j, -i)

		! 	end do
		! end do


		! end function



end module