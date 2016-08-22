module func_var

	use parallel_cubic, Only: parallel
	use interpolation, Only: interp
	use grid_var, Only: g_var

implicit none

	Private
	Public :: f_var

	Type f_var

		Real(8), Allocatable :: h_height(:, :, :)
		Real(8), Allocatable :: x_vel(:, :, :)
		Real(8), Allocatable :: y_vel(:, :, :)
		Real(8), Allocatable :: x_vel_msg(:, :, :)
		Real(8), Allocatable :: y_vel_msg(:, :, :)
		! Real(8), Allocatable :: distance_grid(:, :, :, :)
		real(8) height
		integer(4) step, dim, Xsize, Ysize, interp_factor(1:4), Neighbours_face(6, 4)
		integer(4) ns_x, ns_y, nf_x, nf_y, first_x, first_y, last_x, last_y

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Private :: alloc => alloc
		Procedure, Public :: deinit => deinit
		Procedure, Public :: equal => equal
		Procedure, Public :: start_conditions => start_conditions
		Procedure, Public :: interpolate => interpolate
		Procedure, Public :: Velocity_from_spherical => Velocity_from_spherical
		Procedure, Public :: Velocity_to_spherical => Velocity_to_spherical
	End Type


CONTAINS



	subroutine init(this, paral, height)

		Class(f_var) :: this
		Class(parallel) :: paral
		real(8), intent(in) :: height
		integer(4) :: i

		this.ns_x = paral.ns_xy(1);  this.ns_y = paral.ns_xy(2)
		this.nf_x = paral.nf_xy(1);  this.nf_y = paral.nf_xy(2)

		this.first_x = paral.first_x;  this.first_y = paral.first_y
		this.last_x = paral.last_x;  this.last_y = paral.last_y

		this.Xsize = paral.Xsize;  this.Ysize = paral.Ysize
		this.step = paral.step;  this.height = height;  this.dim = paral.dim

		this.Neighbours_face = paral.Neighbours_face

		this.interp_factor(:) = 0

		do i = 1, 4
			if(paral.Neighbours_face(2, i) /= 2) this.interp_factor(i) = 1
		end do

		call this.alloc()

	end subroutine



	subroutine alloc(this)

		Class(f_var) :: this

		Allocate(this.h_height(this.first_x:this.last_x, this.first_y:this.last_y, 6))
		Allocate(this.x_vel(this.first_x:this.last_x, this.first_y:this.last_y, 6))
		Allocate(this.y_vel(this.first_x:this.last_x, this.first_y:this.last_y, 6))
		Allocate(this.x_vel_msg(this.first_x:this.last_x, this.first_y:this.last_y, 6))
		Allocate(this.y_vel_msg(this.first_x:this.last_x, this.first_y:this.last_y, 6))

	end subroutine



	subroutine deinit(this)
		Class(f_var) :: this

		if (Allocated(this.h_height)) Deallocate(this.h_height)
		if (Allocated(this.x_vel)) Deallocate(this.x_vel)
		if (Allocated(this.y_vel)) Deallocate(this.y_vel)
		if (Allocated(this.x_vel_msg)) Deallocate(this.x_vel_msg)
		if (Allocated(this.y_vel_msg)) Deallocate(this.y_vel_msg)

	end subroutine



	subroutine equal(var_pr, var, grid)

		Class(f_var) :: var_pr, var
		Class(g_var) :: grid

				var_pr.h_height(:, :, :)=var.h_height(:, :, :)
				var_pr.x_vel(:, :, :)=var.x_vel(:, :, :)
				var_pr.y_vel(:, :, :)=var.y_vel(:, :, :)

				var_pr.x_vel_msg(:, :, :)=var_pr.x_vel(:, :, :)
				var_pr.y_vel_msg(:, :, :)=var_pr.y_vel(:, :, :)
				! call var_pr.Velocity_to_spherical(grid)


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
					 h0*exp(-((((10.0/dim)*((x-dim - 0.5)*0.5))**2)+(((10.0/dim)*((y-dim - 0.5)*0.5))**2)))
				end do
			end do
			end if

		end do

	end subroutine




	subroutine interpolate(this, i, g)
		Class(f_var) :: this
		Class(interp) :: i
		Class(g_var) :: g
		! call this.Velocity_from_spherical(g)
		call i.Lagrange(this.h_height, this.interp_factor)
		call i.Lagrange(this.x_vel, this.interp_factor)
		call i.Lagrange(this.y_vel, this.interp_factor)
	end subroutine




	subroutine Velocity_to_spherical(this, g)
		Class(f_var) :: this
		Class(g_var) :: g
		Real(8) :: vel_x_contr, vel_y_contr
		Integer(4) :: x, y, face, i, x_start(4), y_start(4), x_fin(4), y_fin(4)

		x_start(1) = this.ns_x;  x_start(2) = this.nf_x - this.step + 1;  x_start(3) = this.ns_x;  x_start(4) = this.ns_x
		x_fin(1) = this.nf_x;  x_fin(2) = this.nf_x;  x_fin(3) = this.nf_x;  x_fin(4) = this.ns_x + this.step

		y_start(1) = this.nf_y - this.step + 1;  y_start(2) = this.ns_y;  y_start(3) = this.ns_y;  y_start(4) = this.ns_y
		y_fin(1) = this.nf_y;  y_fin(2) = this.nf_y;  y_fin(3) = this.ns_y + this.step;  y_fin(4) = this.nf_y

		do face = 1, 6
			do i = 1, 4

				if(this.interp_factor(i) == 1) then 
					do x = x_start(i), x_fin(i)
						do y = y_start(i), y_fin(i)

							vel_x_contr = g.G_inverse(x, y, 1, 1) * this.x_vel(x, y, face) + g.G_inverse(x, y, 1, 2) * this.y_vel(x, y, face)
							vel_y_contr = g.G_inverse(x, y, 2, 2) * this.y_vel(x, y, face) + g.G_inverse(x, y, 2, 1) * this.x_vel(x, y, face)

							this.x_vel_msg(x, y, face) = g.To_sph_coord(1, 1, x, y, face) * vel_x_contr + g.To_sph_coord(1, 2, x, y, face) * vel_y_contr
							this.y_vel_msg(x, y, face) = g.To_sph_coord(2, 1, x, y, face) * vel_x_contr + g.To_sph_coord(2, 2, x, y, face) * vel_y_contr

						end do
					end do
				end if

			end do
		end do

	end subroutine



	subroutine Velocity_from_spherical(this, g)
		Class(f_var) :: this
		Class(g_var) :: g
		Real(8) :: vel_x_contr, vel_y_contr
		Integer(4) :: x, y, face, i, neib_face, x_start(4), y_start(4), x_fin(4), y_fin(4)

		x_start(1) = this.ns_x;  x_start(2) = this.nf_x + 1;  x_start(3) = this.ns_x;  x_start(4) = this.ns_x - this.step
		x_fin(1) = this.nf_x;  x_fin(2) = this.nf_x + this.step;  x_fin(3) = this.nf_x;  x_fin(4) = this.ns_x - 1

		y_start(1) = this.nf_y + 1;  y_start(2) = this.ns_y;  y_start(3) = this.ns_y - this.step;  y_start(4) = this.ns_y
		y_fin(1) = this.nf_y + this.step;  y_fin(2) = this.nf_y;  y_fin(3) = this.ns_y - 1;  y_fin(4) = this.nf_y

		do face = 1, 6
			do i = 1, 4

				if(this.interp_factor(i) == 1) then 
					do x = x_start(i), x_fin(i)
						do y = y_start(i), y_fin(i)

							vel_x_contr = g.From_sph_coord(1, 1, x, y, face) * this.x_vel_msg(x, y, face) + g.From_sph_coord(1, 2, x, y, face) * this.y_vel_msg(x, y, face)
							vel_y_contr = g.From_sph_coord(2, 1, x, y, face) * this.x_vel_msg(x, y, face) + g.From_sph_coord(2, 2, x, y, face) * this.y_vel_msg(x, y, face)

							this.x_vel(x, y, face) = g.G_tensor(x, y, 1, 1) * vel_x_contr + g.G_tensor(x, y, 1, 2) * vel_y_contr
							this.y_vel(x, y, face) = g.G_tensor(x, y, 2, 1) * vel_x_contr + g.G_tensor(x, y, 2, 2) * vel_y_contr

						end do
					end do
				end if

			end do
		end do

	end subroutine




end module