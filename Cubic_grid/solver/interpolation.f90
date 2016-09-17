module interpolation

	use grid_var, Only: g_var
	use sphere_geometry, Only: geometry
	use mpi

	implicit none


	Private
	Public :: interp

	Type interp

		Real(8), Allocatable :: weight(:, :, :)
		Integer(4), Allocatable :: x0_mass(:, :)
		integer(4) ns_x, ns_y, nf_x, nf_y, n, dim, step, first_x, first_y, last_x, last_y, snd_xy(6, 4, 2), rcv_xy(6, 4, 2)

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Public :: Lagrange => Lagrange
		Procedure, Private :: X0_find => X0_find
		Procedure, Private :: weight_find => weight_find
	End Type

	CONTAINS


	Subroutine init(this, g, n)

		Class(interp) :: this
		Class(g_var) :: g
		Integer(4) :: n
		Integer(4) :: k, i, j, xk, x, x0, dim

		this.ns_x = g.ns_xy(1);  this.nf_x = g.nf_xy(1);  this.ns_y = g.ns_xy(2);  this.nf_y = g.nf_xy(2)
		this.n = n;  this.dim = g.dim;  dim = this.dim; this.step = g.step

		this.first_x = g.first_x;  this.first_y = g.first_y
		this.last_x = g.last_x;  this.last_y = g.last_y

		this.snd_xy(:,:,:) = g.snd_xy(:,:,:)
		this.rcv_xy(:,:,:) = g.rcv_xy(:,:,:)

		Allocate(this.x0_mass(1: 2*dim, 1: g.step))
		Allocate(this.weight(1: 2*dim, 1: 2*dim, 1: g.step))

		do i = 1, g.step
			do x = 1, 2*dim
				call this.X0_find(g, x, i, x0)
				this.x0_mass(x, i) = x0
			end do
		end do

		do i = 1, g.step
			do x = 1, 2*dim
				do xk = this.x0_mass(x, i) - n/2, this.x0_mass(x, i) + (n+1)/2
					this.weight(xk, x, i) = this.weight_find(g, xk, x, i)
				end do
			end do
		end do

	End Subroutine



	Subroutine Lagrange(this, Mass, interp_factor)

		Class(interp) :: this

		Real(8), Intent(inout) :: Mass(this.first_x:this.last_x, this.first_y:this.last_y, 6)
		Integer(4), Intent(in) :: interp_factor(1:4)
		Real(8) :: temp(1:6), Mass_temp(this.first_x:this.last_x, this.first_y:this.last_y, 6)
		Integer(4) xj, xk, k, j, n, x0, set_x, x, s, y, yk, i, face, x_fin(4), y_fin(4)

		x_fin(:) = this.last_x;  x_fin(4) = this.rcv_xy(2, 4, 1) + this.step
		y_fin(:) = this.last_y;  y_fin(3) = this.rcv_xy(2, 3, 2) + this.step

		n=this.n
		Mass_temp(:, :, :) = Mass(:, :, :)

		do face = 1, 6

		if(interp_factor(1) == 1) then !up
			do x = this.ns_x, this.nf_x
				do y = this.nf_y + 1, this.nf_y + this.step
					x0 = this.x0_mass(x, y - this.nf_y)
					s = x0 - n/2
					temp=0.0

					do xk=s, s + n
						temp(face) = temp(face) + (Mass(xk, y, face)*this.weight(xk, x, y - this.nf_y))
					end do

					Mass_temp(x, y, face) = temp(face)
				end do
			end do
		end if

			Mass(:, :, :) = Mass_temp(:, :, :)
			Mass_temp(:, :, :) = Mass(:, :, :)

		if(interp_factor(2) == 1) then !right
			do y = this.ns_y, this.nf_y
				do x = this.nf_x + 1, this.nf_x + this.step
					x0 = this.x0_mass(y, x - this.nf_x)
					s = x0 - n/2
					temp=0.0

					do yk=s, s + n
						temp(face) = temp(face) + (Mass(x, yk, face)*this.weight(yk, y, x - this.nf_x))
					end do

					Mass_temp(x, y, face) = temp(face)
				end do
			end do
		end if

		Mass(:, :, :) = Mass_temp(:, :, :)
		Mass_temp(:, :, :) = Mass(:, :, :)

		if(interp_factor(3) == 1) then !down
			do x = this.ns_x, this.nf_x
				do y = this.ns_y - 1, this.ns_y - this.step
					x0 = this.x0_mass(x, this.ns_y - y)
					s = x0 - n/2
					temp=0.0

					do xk=s, s + n
						temp(face) = temp(face) + (Mass(xk, y, face)*this.weight(xk, x, this.ns_y - y))
					end do

					Mass_temp(x, y, face) = temp(face)
				end do
			end do
		end if

		Mass(:, :, :) = Mass_temp(:, :, :)
		Mass_temp(:, :, :) = Mass(:, :, :)

		if(interp_factor(4) == 1) then !left
			do y = this.ns_y, this.nf_y
				do x = this.ns_x - 1, this.ns_x - this.step
					x0 = this.x0_mass(y, this.ns_x - x)
					s = x0 - n/2
					temp=0.0

					do yk=s, s + n
						temp(face) = temp(face) + (Mass(x, yk, face)*this.weight(yk, y, this.ns_x - x))
					end do

					Mass_temp(x, y, face) = temp(face)
				end do
			end do
		end if

		end do

		Mass(:, :, :) = Mass_temp(:, :, :)

	End Subroutine



	Subroutine X0_find(this, g, x, step, x0)

		Class(interp) :: this
		Class(g_var) :: g
		Integer(4), intent(in) :: x, step
		Integer(4), intent(out) :: x0
		Integer(4) :: xk
		Real(8) :: gap

		do xk = 1, 2*this.dim

			gap = (g.latlon_c(1, 1 - step, x, 2)) - (g.latlon_c(1, step, xk, 2))
			if(gap >= 0.0) then
				x0 = xk
			end if

		end do

	end Subroutine



	Real(8) function weight_find(this, g, xk, x, step)

		Class(interp) :: this
		Class(g_var) :: g
		Type(geometry) :: geom
		Integer(4), intent(in) :: step, x, xk
		Integer(4) :: xj, x0, gap, n, s
		Real(8) :: temp, latlon_x(1:2), latlon_xk(1:2), latlon_xj(1:2)

		n=this.n

		x0 = this.x0_mass(x, step)
		temp=1.0;  s = x0 - n/2

			do xj=s, s + (n+1)/2
				if(xj /= xk)then
					latlon_x(1) = g.latlon_c(1, 1-step, x, 2)
					latlon_x(2) = g.latlon_c(2, step, xk, 2)
					latlon_xk(:) = g.latlon_c(:, step, xk, 2)
					latlon_xj(:) = g.latlon_c(:, step, xj, 2)
					temp = (geom.angle(latlon_x(:), latlon_xk(:)))/(geom.angle(latlon_xj(:), latlon_xk(:)))
! 					temp = (g.cube_coord_c(2, 1-step, x) - g.cube_coord_c(2, step, xk))/(g.cube_coord_c(2, step, xk) - g.cube_coord_c(2, step, xj))
				end if
			end do

		weight_find = temp
		! print *, temp, x


	end function



end module
