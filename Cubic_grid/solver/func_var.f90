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
	End Type


CONTAINS



	subroutine init(this, paral, height)

		Class(f_var) :: this
		Class(parallel) :: paral
		real(8), intent(in) :: height

		this.ns_x = paral.ns_xy(1);  this.ns_y = paral.ns_xy(2)
		this.nf_x = paral.nf_xy(1);  this.nf_y = paral.nf_xy(2)

		this.first_x = paral.first_x;  this.first_y = paral.first_y
		this.last_x = paral.last_x;  this.last_y = paral.last_y

		this.Xsize = paral.Xsize;  this.Ysize = paral.Ysize
		this.step = paral.step;  this.height = height;  this.dim = paral.dim

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

				var_pr.h_height(:, :, :)=var.h_height(:, :, :)
				var_pr.x_vel(:, :, :)=var.x_vel(:, :, :)
				var_pr.y_vel(:, :, :)=var.y_vel(:, :, :)

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

			if ( face == 4 ) then
			do y = this.first_y, this.last_y
				do x = this.first_x, this.last_x
					this.h_height(x, y, face) =&
					 h0*exp(-((((10.0/dim)*((x-dim)*0.5))**2)+(((10.0/dim)*((y-dim)*0.5))**2)))
				end do
			end do
			end if

		end do

	end subroutine



end module