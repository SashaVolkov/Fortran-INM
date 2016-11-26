module interpolation

	use grid_var, Only: g_var
	use sphere_geometry, Only: geometry
	use mpi
	use omp_lib

	implicit none


	Private
	Public :: interp

	Type interp

		Real(8), Allocatable :: weight(:, :, :, :)
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

		Allocate(this.x0_mass(2*dim, g.step))
		Allocate(this.weight(n, 2*dim, 2*dim, g.step))

		do i = 1, g.step
			do x = 1, 2*dim
				call this.X0_find(g, x, i, x0)
				this.x0_mass(x, i) = x0
			end do
		end do

		do i = 1, g.step
			do x = 1, 2*dim
				x0 = this.x0_mass(x, i)
				call this.weight_find(g, x0, x, i, this.weight(:, x0, x, i))
				! print *, i, x, x0, real(this.weight(:, x0, x, i),4)
			end do
		end do

	End Subroutine



	Subroutine Lagrange(this, Mass, interp_factor)

		Class(interp) :: this

		Real(8), Intent(inout) :: Mass(this.first_x:this.last_x, this.first_y:this.last_y, 6)
		Integer(4), Intent(in) :: interp_factor(1:4)
		Real(8) :: Mass_temp(this.first_x:this.last_x, this.first_y:this.last_y, 6)
		Integer(4) k, n, x0, x, y, i, face, x_fin(4), y_fin(4), x_int(1:this.n), dim

		x_fin(:) = this.last_x;  x_fin(4) = this.rcv_xy(2, 4, 1) + this.step
		y_fin(:) = this.last_y;  y_fin(3) = this.rcv_xy(2, 3, 2) + this.step

		n=this.n;  dim = this.dim
		Mass_temp(:, :, :) = Mass(:, :, :)

		do face = 1, 6

		if(interp_factor(1) == 1) then !up
		y = 2*this.dim
			do x = this.ns_x, this.nf_x
				do i = 1, this.step

					x0 = this.x0_mass(x, i);  n=this.n
					if(i == 1 .and. (x==this.ns_x .or. x==this.nf_x)) n = 2
					x_int(1) = x0 - n/2 + 1
					do k = 2, n
						x_int(k) = x_int(k-1)+1
					end do

					Mass_temp(x, y+i, face) = 0d0
					do k = 1, n
						Mass_temp(x, y+i, face) = Mass(x_int(k), y+i, face)*this.weight(k, x0, x, i) + Mass_temp(x, y+i, face)
					end do

				end do
			end do
			Mass(:, :, face) = Mass_temp(:, :, face)
			Mass_temp(:, :, face) = Mass(:, :, face)
		end if


		if(interp_factor(2) == 1) then !right
		x = 2*this.dim
			do y = this.ns_y, this.nf_y
				do i = 1, this.step

					x0 = this.x0_mass(y, i);  n=this.n
					if(i == 1 .and. (y==this.ns_y .or. y==this.nf_y)) n = 2
					x_int(1) = x0 - n/2 + 1
					do k = 2, n
						x_int(k) = x_int(k-1)+1
					end do
					Mass_temp(x+i, y, face) = 0d0
					do k = 1, n
						Mass_temp(x+i, y, face) = Mass(x+i, x_int(k), face)*this.weight(k, x0, y, i) + Mass_temp(x+i, y, face)
					end do

				end do
			end do
			Mass(:, :, face) = Mass_temp(:, :, face)
			Mass_temp(:, :, face) = Mass(:, :, face)
		end if


		if(interp_factor(3) == 1) then !down
			y = 1
			do x = this.ns_x, this.nf_x
				do i = 1, this.step

					x0 = this.x0_mass(x, i);  n=this.n
					if(i == 1 .and. (x==this.ns_x .or. x==this.nf_x)) n = 2
					x_int(1) = x0 - n/2 + 1
					do k = 2, n
						x_int(k) = x_int(k-1)+1
					end do

					Mass_temp(x, y-i, face) = 0d0
					do k = 1, n
						Mass_temp(x, y-i, face) = Mass(x_int(k), y-i, face)*this.weight(k, x0, x, i) + Mass_temp(x, y-i, face)
					end do
				end do
			end do
			Mass(:, :, face) = Mass_temp(:, :, face)
			Mass_temp(:, :, face) = Mass(:, :, face)
		end if


		if(interp_factor(4) == 1) then !left
		x = 1
			do y = this.ns_y, this.nf_y
				do i = 1, this.step

					x0 = this.x0_mass(y, i);  n=this.n
					if(i == 1 .and. (y==this.ns_y .or. y==this.nf_y)) n = 2
					x_int(1) = x0 - n/2 + 1
					do k = 2, n
						x_int(k) = x_int(k-1)+1
					end do
					Mass_temp(x-i, y, face) = 0d0
					do k = 1, n
						Mass_temp(x-i, y, face) = Mass(x-i, x_int(k), face)*this.weight(k, x0, y, i) + Mass_temp(x-i, y, face)
					end do
				end do
			end do
			Mass(:, :, face) = Mass_temp(:, :, face)
			Mass_temp(:, :, face) = Mass(:, :, face)
		end if

		end do


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



	Subroutine weight_find(this, g, x0, x, step, weight)

		Class(interp) :: this
		Class(g_var) :: g
		Type(geometry) :: geom
		Integer(4), intent(in) :: step, x, x0
		Real(8), intent(out) :: weight(1:this.n)
		Integer(4) :: n, y(1:this.n), i, j, dim
		Real(8) :: s, numen, k(4, 4), h(1:4), latlon_y(this.n), latlon_x

		n = this.n;  dim = this.dim; weight(:) = 1d0
		if(step == 1 .and. (x==this.ns_y .or. x==this.nf_y .or. x==this.ns_x .or. x==this.nf_x)) n = 2
		latlon_x = g.latlon_c(1,1-step,x,2)

		y(1) = x0 - n/2 + 1
		latlon_y(1) = g.latlon_c(1,2*g.dim+1-step,y(1),5)
		h(1) = abs(latlon_x - latlon_y(1))

		do i = 2, n
			y(i) = y(i-1)+1
			latlon_y(i) = g.latlon_c(1,2*g.dim+1-step,y(i),5)
			h(i) = latlon_x - latlon_y(i)
			do j = 1, n
				if(j /= i) k(i, j) = h(i)/(latlon_y(j) - latlon_y(i))
				if(j /= i) k(j, i) = h(j)/(latlon_y(i) - latlon_y(j))
			end do
		end do

		do i = 1, n
			do j = 1, n
				if(j /= i) weight(i) = weight(i)*k(j,i)
			end do
		end do

	end Subroutine



end module
