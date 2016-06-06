module func_var

	use parallel_cubic, Only: parallel

implicit none

	Private
	Public :: f_var

	Type f_var

		Real(8), Allocatable :: h_height(:, :)
		Real(8), Allocatable :: u_vel(:, :)
		Real(8), Allocatable :: v_vel(:, :)
		! Real(8), Allocatable :: distance_grid(:, :, :, :)
		real(8) height
		integer(4) step, face, dim, Xsize, Ysize
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



	subroutine init(this, paral, dim, step, height, face)

		Class(f_var) :: this
		Class(parallel) :: paral
		integer(4), intent(in) :: step, face, dim
		real(8), intent(in) :: height


		this.ns_x = paral.ns_xy(1);  this.ns_y = paral.ns_xy(2)
		this.nf_x = paral.nf_xy(1);  this.nf_y = paral.nf_xy(2)

		this.first_x = paral.ns_xy(1) - step;  this.first_y = paral.ns_xy(2) - step
		this.last_x = paral.nf_xy(1) + step;  this.last_y = paral.nf_xy(2) + step

		this.Xsize = paral.Xsize;  this.Ysize = paral.Ysize

		this.step = step;  this.height = height;  this.face = face; this.dim = dim

		call this.alloc()

	end subroutine



	subroutine alloc(this)

		Class(f_var) :: this

		Allocate(this.h_height(this.first_x:this.last_x, this.first_y:this.last_y))
		Allocate(this.u_vel(this.first_x:this.last_x, this.first_y:this.last_y))
		Allocate(this.v_vel(this.first_x:this.last_x, this.first_y:this.last_y))

	end subroutine



	subroutine deinit(this)
		Class(f_var) :: this

		if (Allocated(this.h_height)) Deallocate(this.h_height)
		if (Allocated(this.u_vel)) Deallocate(this.u_vel)
		if (Allocated(this.v_vel)) Deallocate(this.v_vel)

	end subroutine



	subroutine equal(var_pr, var)

		Class(f_var) :: var_pr, var
		integer(4) dim, x, y

			! do y = var_pr.first_y, var_pr.last_y
			! 	do x = var_pr.first_x, var_pr.last_x
				var_pr.h_height(:, :)=var.h_height(:, :)
				var_pr.u_vel(:, :)=var.u_vel(:, :)
				var_pr.v_vel(:, :)=var.v_vel(:, :)
			! 	end do
			! end do


	end subroutine



	subroutine start_conditions(this)

		Class(f_var) :: this
		integer(4) dim
		real(8) h0

		integer(4) x, y

		h0 = this.height;  dim = this.dim


		do y = this.first_y, this.last_y
			do x = this.first_x, this.last_x
				this.h_height(x, y) = 0
				this.u_vel(x, y) = 0
				this.v_vel(x, y) = 0
			end do
		end do

		if ( this.face == 2 ) then
		do y = this.first_y, this.last_y
			do x = this.first_x, this.last_x
				this.h_height(x, y) =&
				 h0*exp(-((((10.0/dim)*(x*0.5))**2)+(((10.0/dim)*(y*0.5))**2)))
			end do
		end do
		end if


	end subroutine



		! integer(4) function borders(this, var, var_grey)

		! Class(f_var) :: this
		! Class(f_var) :: var(1:6), var_grey(1:6)
		! integer(4) dim, step, dim_st, x, y, i, j

		! dim = var(1).dim;  step = var(1).step;  dim_st = var(1).dim_st


		! do y = 0, step
		! 	do x = -dim, dim

		! 		var_grey(1).h_height(x, -dim - y) = var(2).h_height(x, dim - y)
		! 		var_grey(1).u_vel(x, -dim - y) = var(2).u_vel(x, dim - y)
		! 		var_grey(1).v_vel(x, -dim - y) = var(2).v_vel(x, dim - y)

		! 		var_grey(2).h_height(x, dim + y) = var(1).h_height(x, -dim + y)
		! 		var_grey(2).u_vel(x, dim + y) = var(1).u_vel(x, -dim + y)
		! 		var_grey(2).v_vel(x, dim + y) = var(1).v_vel(x, -dim + y)


		! 		var_grey(6).h_height(x, dim + y) = var(2).h_height(x, -dim + y)
		! 		var_grey(6).u_vel(x, dim + y) = var(2).u_vel(x, -dim + y)
		! 		var_grey(6).v_vel(x, dim + y) = var(2).v_vel(x, -dim + y)

		! 		var_grey(2).h_height(x, -dim - y) = var(6).h_height(x, dim - y)
		! 		var_grey(2).u_vel(x, -dim - y) = var(6).u_vel(x, dim - y)
		! 		var_grey(2).v_vel(x, -dim - y) = var(6).v_vel(x, dim - y)

		! 	end do
		! end do


		! do y = -dim, dim
		! 	do x = 0, step

		! 		var_grey(5).h_height(dim + x, y) = var(2).h_height(-dim + x, y)
		! 		var_grey(5).u_vel(dim + x, y) = var(2).u_vel(-dim + x, y)
		! 		var_grey(5).v_vel(dim + x, y) = var(2).v_vel(-dim + x, y)

		! 		var_grey(2).h_height(-dim - x, y) = var(5).h_height(dim - x, y)
		! 		var_grey(2).u_vel(-dim - x, y) = var(5).u_vel(dim - x, y)
		! 		var_grey(2).v_vel(-dim - x, y) = var(5).v_vel(dim - x, y)


		! 		var_grey(3).h_height(-dim - x, y) = var(2).h_height(dim - x, y)
		! 		var_grey(3).u_vel(-dim - x, y) = var(2).u_vel(dim - x, y)
		! 		var_grey(3).v_vel(-dim - x, y) = var(2).v_vel(dim - x, y)

		! 		var_grey(2).h_height(dim + x, y) = var(3).h_height(-dim + x, y)
		! 		var_grey(2).u_vel(dim + x, y) = var(3).u_vel(-dim + x, y)
		! 		var_grey(2).v_vel(dim + x, y) = var(3).v_vel(-dim + x, y)

		! 	end do
		! end do


		! do y = 0, step
		! 	do x = -dim, dim

		! 	var_grey(1).h_height(-x, dim + y) = var(4).h_height(x, dim - y)
		! 	var_grey(1).u_vel(-x, dim + y) = -var(4).u_vel(x, dim - y)
		! 	var_grey(1).v_vel(-x, dim + y) = -var(4).v_vel(x, dim - y)

		! 		var_grey(4).h_height(x, dim + y) = var(1).h_height(-x, dim - y)
		! 		var_grey(4).u_vel(x, dim + y) = -var(1).u_vel(-x, dim - y)
		! 		var_grey(4).v_vel(x, dim + y) = -var(1).v_vel(-x, dim - y)


		! 		var_grey(6).h_height(-x, -dim - y) = var(4).h_height(x, -dim + y)
		! 		var_grey(6).u_vel(-x, -dim - y) = -var(4).u_vel(x, -dim + y)
		! 		var_grey(6).v_vel(-x, -dim - y) = -var(4).v_vel(x, -dim + y)

		! 		var_grey(4).h_height(x, -dim - y) = var(6).h_height(-x, -dim + y)
		! 		var_grey(4).u_vel(x, -dim - y) = -var(6).u_vel(-x, -dim + y)
		! 		var_grey(4).v_vel(x, -dim - y) = -var(6).v_vel(-x, -dim + y)

		! 	end do
		! end do


		! do y = -dim, dim
		! 	do x = 0, step

		! 		var_grey(3).h_height(dim + x, y) = var(4).h_height(-dim + x, y)
		! 		var_grey(3).u_vel(dim + x, y) = var(4).u_vel(-dim + x, y)
		! 		var_grey(3).v_vel(dim + x, y) = var(4).v_vel(-dim + x, y)

		! 		var_grey(4).h_height(-dim - x, y) = var(3).h_height(dim - x, y)
		! 		var_grey(4).u_vel(-dim - x, y) = var(3).u_vel(dim - x, y)
		! 		var_grey(4).v_vel(-dim - x, y) = var(3).v_vel(dim - x, y)


		! 		var_grey(4).h_height(dim + x, y) = var(5).h_height(-dim + x, y)
		! 		var_grey(4).u_vel(dim + x, y) = var(5).u_vel(-dim + x, y)
		! 		var_grey(4).v_vel(dim + x, y) = var(5).v_vel(-dim + x, y)

		! 		var_grey(5).h_height(-dim - x, y) = var(4).h_height(dim - x, y)
		! 		var_grey(5).u_vel(-dim - x, y) = var(4).u_vel(dim - x, y)
		! 		var_grey(5).v_vel(-dim - x, y) = var(4).v_vel(dim - x, y)

		! 	end do
		! end do


		! do i = -dim, dim
		! 	do j = 0, step

		! 		var_grey(1).h_height(-dim - j, -i) = var(5).h_height(i, dim - j)
		! 		var_grey(1).u_vel(-dim - j, -i) = var(5).v_vel(i, dim - j)
		! 		var_grey(1).v_vel(-dim - j, -i) = -var(5).u_vel(i, dim - j)

		! 		var_grey(5).h_height(i, dim + j) = var(1).h_height(-dim + j, -i)
		! 		var_grey(5).u_vel(i, dim + j) = -var(1).v_vel(-dim + j, -i)
		! 		var_grey(5).v_vel(i, dim + j) = var(1).u_vel(-dim + j, -i)


		! 		var_grey(6).h_height(-dim - j, i) = var(5).h_height(i, -dim + j)
		! 		var_grey(6).u_vel(-dim - j, i) = -var(5).v_vel(i, -dim + j)
		! 		var_grey(6).v_vel(-dim - j, i) = var(5).u_vel(i, -dim + j)

		! 		var_grey(5).h_height(i, -dim - j) = var(6).h_height(-dim + j, i)
		! 		var_grey(5).u_vel(i, -dim - j) = var(6).v_vel(-dim + j, i)
		! 		var_grey(5).v_vel(i, -dim - j) = -var(6).u_vel(-dim + j, i)



		! 		var_grey(1).h_height(dim + j, i) = var(3).h_height(i, dim - j)
		! 		var_grey(1).u_vel(dim + j, i) = -var(3).v_vel(i, dim - j)
		! 		var_grey(1).v_vel(dim + j, i) = var(3).u_vel(i, dim - j)

		! 		var_grey(3).h_height(i, dim + j) = var(1).h_height(dim - j, i)
		! 		var_grey(3).u_vel(i, dim + j) = var(1).v_vel(dim - j, i)
		! 		var_grey(3).v_vel(i, dim + j) = -var(1).u_vel(dim - j, i)


		! 		var_grey(6).h_height(dim + j, -i) = var(3).h_height(i, -dim + j)
		! 		var_grey(6).u_vel(dim + j, -i) = var(3).v_vel(i, -dim + j)
		! 		var_grey(6).v_vel(dim + j, -i) = -var(3).u_vel(i, -dim + j)

		! 		var_grey(3).h_height(i, -dim - j) = var(6).h_height(dim - j, -i)
		! 		var_grey(3).u_vel(i, -dim - j) = -var(6).v_vel(dim - j, -i)
		! 		var_grey(3).v_vel(i, -dim - j) = var(6).u_vel(dim - j, -i)

		! 	end do
		! end do


		! end function



end module