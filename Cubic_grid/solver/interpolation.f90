module interpolation

	use metrics, Only: metric
	use sphere_geometry, Only: geometry
	use mpi
	use omp_lib

	implicit none


	Private
	Public :: interp

	Type interp

		Real(8), Allocatable :: weight(:, :, :, :)
		Real(8) :: fantom1(2)
		Integer(4), Allocatable :: x0_mass(:, :)
		Integer(4) ns_x, ns_y, nf_x, nf_y, n, dim, step, first_x, first_y, last_x, last_y, snd_xy(6, 4, 2), rcv_xy(6, 4, 2)

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Public :: Lagrange => Lagrange
		Procedure, Private :: X0_find => X0_find
		Procedure, Private :: weight_find => weight_find
	End Type

	CONTAINS


	Subroutine init(this, metr, n)

		Class(interp) :: this
		Class(metric) :: metr
		Integer(4) :: n
		Integer(4) :: k, i, j, xk, x, x0, dim, step

		this.ns_x = metr.ns_xy(1);  this.nf_x = metr.nf_xy(1);  this.ns_y = metr.ns_xy(2);  this.nf_y = metr.nf_xy(2)
		this.n = n;  this.dim = metr.dim;  dim = this.dim;  this.step = metr.step;  step = this.step

		this.first_x = metr.first_x;  this.first_y = metr.first_y
		this.last_x = metr.last_x;  this.last_y = metr.last_y

		this.snd_xy = metr.snd_xy;  this.rcv_xy = metr.rcv_xy

		Allocate(this.x0_mass(2*dim, step))
		Allocate(this.weight(n, 2*dim, 2*dim, step))

		do i = 1, step
			do x = 1, 2*dim
				call this.X0_find(metr, x, i, x0)
				this.x0_mass(x, i) = x0
			end do
		end do

		do i = 1, step
			do x = 1, 2*dim
				x0 = this.x0_mass(x, i)
				call this.weight_find(metr, x0, x, i, this.weight(:, x0, x, i))
			end do
		end do

	End Subroutine



	Subroutine Lagrange(this, Mass, interp_factor)

		Class(interp) :: this

		Real(8), Intent(inout) :: Mass(this.first_x:this.last_x, this.first_y:this.last_y, 6)
		Integer(4), Intent(in) :: interp_factor(1:4)
		Real(8) :: Mass_temp(this.first_x:this.last_x, this.first_y:this.last_y, 6), M, M1, M2
		Integer(4) k, n, x0, x, y, i, face, x_fin(4), y_fin(4), x_int(1:this.n), dim

		x_fin(:) = this.last_x;  x_fin(4) = this.rcv_xy(2, 4, 1) + this.step
		y_fin(:) = this.last_y;  y_fin(3) = this.rcv_xy(2, 3, 2) + this.step

		n=this.n;  dim = this.dim
		Mass_temp = Mass

		do face = 1, 6

		if(interp_factor(1) == 1) then !up
		y = 2*this.dim
			do x = this.ns_x, this.nf_x
				do i = 1, this.step

					x0 = this.x0_mass(x, i);  n=this.n
					x_int(1) = x0 - n/2 + 1
					do k = 2, n
						x_int(k) = x_int(k-1)+1
					end do

					M2 = 0d0
					do k = 1, n
						M = Mass(x_int(k), y+i, face)

						if(x_int(k) < 1 .and. x0 == 1 .and. i==1) then
							! M1 = Mass(2*dim, 1, face)*this.fantom1(1) + Mass(2*dim-1, 1, face)*this.fantom1(2)
							! M = M1*this.weight(1,1,1,1) + Mass(0,2*dim,face)*this.weight(2,1,1,1) + Mass(-1,2*dim,face)*this.weight(3,1,1,1) + Mass(-2,2*dim,face)*this.weight(4,1,1,1)
							M = Mass(0,2*dim,face)*this.fantom1(1) + Mass(-1,2*dim,face)*this.fantom1(2)
						end if

						if(x_int(k) > 2*dim .and. x0 == 2*dim-1 .and. i==1) then
							! M1 = Mass(2*dim, 2*dim-1, face)*this.fantom1(1) + Mass(2*dim, 2*dim, face)*this.fantom1(2)
							! M = M1*this.weight(1,1,1,1) + Mass(2*dim+1,2*dim,face)*this.weight(2,1,1,1) + Mass(2*dim+2,2*dim,face)*this.weight(3,1,1,1) + Mass(2*dim+3,2*dim,face)*this.weight(4,1,1,1)
							M = Mass(2*dim+1,2*dim,face)*this.fantom1(1) + Mass(2*dim+2,2*dim,face)*this.fantom1(2)
						end if

						M2 = M*this.weight(k, x0, x, i) + M2
					end do
					Mass_temp(x, y+i, face) = M2
				end do
			end do
		end if


		if(interp_factor(2) == 1) then !right
		x = 2*this.dim
			do y = this.ns_y, this.nf_y
				do i = 1, this.step

					x0 = this.x0_mass(y, i);  n=this.n
					x_int(1) = x0 - n/2 + 1
					do k = 2, n
						x_int(k) = x_int(k-1)+1
					end do
					M2 = 0d0
					do k = 1, n
						M = Mass(x+i, x_int(k), face)

						if(x_int(k) < 1 .and. x0 == 1 .and. i==1) then
							! M1 = Mass(2*dim, 1, face)*this.fantom1(1) + Mass(2*dim-1, 1, face)*this.fantom1(2)
							! M = M1*this.weight(1,1,1,1) + Mass(2*dim,0,face)*this.weight(2,1,1,1) + Mass(2*dim,-1,face)*this.weight(3,1,1,1) + Mass(2*dim,-2,face)*this.weight(4,1,1,1)
							M = Mass(2*dim,0,face)*this.fantom1(1) + Mass(2*dim,-1,face)*this.fantom1(2)
						end if

						if(x_int(k) > 2*dim .and. x0 == 2*dim-1 .and. i==1) then
							! M1 = Mass(2*dim-1, 2*dim, face)*this.fantom1(1) + Mass(2*dim, 2*dim, face)*this.fantom1(2)
							! M = M1*this.weight(1,1,1,1) + Mass(2*dim,2*dim+1,face)*this.weight(2,1,1,1) + Mass(2*dim,2*dim+2,face)*this.weight(3,1,1,1) + Mass(2*dim,2*dim+3,face)*this.weight(4,1,1,1)
							M = Mass(2*dim,2*dim+1,face)*this.fantom1(1) + Mass(2*dim,2*dim+2,face)*this.fantom1(2)
						end if

						M2 = M*this.weight(k, x0, y, i) + M2
					end do
					Mass_temp(x+i, y, face) = M2
				end do
			end do
		end if


		if(interp_factor(3) == 1) then !down
			y = 1
			do x = this.ns_x, this.nf_x
				do i = 1, this.step

					x0 = this.x0_mass(x, i);  n=this.n
					x_int(1) = x0 - n/2 + 1
					do k = 2, n
						x_int(k) = x_int(k-1)+1
					end do

					M2 = 0d0
					do k = 1, n
						M = Mass(x_int(k), y-i, face)

						if(x_int(k) < 1 .and. x0 == 1 .and. i==1) then
							! M1 = Mass(1, 1, face)*this.fantom1(1) + Mass(1, 2, face)*this.fantom1(2)
							! M = M1*this.weight(1,1,1,1) + Mass(0,1,face)*this.weight(2,1,1,1) + Mass(-1,1,face)*this.weight(3,1,1,1) + Mass(-2,1,face)*this.weight(4,1,1,1)
							M = Mass(0,1,face)*this.fantom1(1) + Mass(-1,1,face)*this.fantom1(2)
						end if

						if(x_int(k) > 2*dim .and. x0 == 2*dim-1 .and. i==1) then
							! M1 = Mass(2*dim, 1, face)*this.fantom1(1) + Mass(2*dim, 2, face)*this.fantom1(2)
							! M = M1*this.weight(1,1,1,1) + Mass(2*dim+1,1,face)*this.weight(2,1,1,1) + Mass(2*dim+2,1,face)*this.weight(3,1,1,1) + Mass(2*dim+3,1,face)*this.weight(4,1,1,1)
							M = Mass(2*dim+1,1,face)*this.fantom1(1) + Mass(2*dim+2,1,face)*this.fantom1(2)
						end if

						M2 = M*this.weight(k, x0, x, i) + M2
					end do
					Mass_temp(x, y-i, face) = M2
				end do
			end do
		end if


		if(interp_factor(4) == 1) then !left
		x = 1
			do y = this.ns_y, this.nf_y
				do i = 1, this.step

					x0 = this.x0_mass(y, i);  n=this.n
					x_int(1) = x0 - n/2 + 1
					do k = 2, n
						x_int(k) = x_int(k-1)+1
					end do
					M2 = 0d0
					do k = 1, n
						M = Mass(x-i, x_int(k), face)

						if(x_int(k) < 1 .and. x0 == 1 .and. i==1) then
							! M1 = Mass(1, 1, face)*this.fantom1(1) + Mass(2, 1, face)*this.fantom1(2)
							! M = M1*this.weight(1,1,1,1) + Mass(1,0,face)*this.weight(2,1,1,1) + Mass(1,-1,face)*this.weight(3,1,1,1) + Mass(1,-2,face)*this.weight(4,1,1,1)
							M = Mass(1,0,face)*this.fantom1(1) + Mass(1,-1,face)*this.fantom1(2)
						end if

						if(x_int(k) > 2*dim .and. x0 == 2*dim-1 .and. i==1) then
							! M1 = Mass(1, 2*dim, face)*this.fantom1(1) + Mass(2, 2*dim, face)*this.fantom1(2)
							! M = M1*this.weight(1,1,1,1) + Mass(1,2*dim+1,face)*this.weight(2,1,1,1) + Mass(1,2*dim+2,face)*this.weight(3,1,1,1) + Mass(1,2*dim+3,face)*this.weight(4,1,1,1)
							M = Mass(1,2*dim+1,face)*this.fantom1(1) + Mass(1,2*dim+2,face)*this.fantom1(2)
						end if

						M2 = M*this.weight(k, x0, y, i) + M2
					end do
					Mass_temp(x-i, y, face) = M2
				end do
			end do
		end if

		end do

		Mass = Mass_temp

	End Subroutine



	Subroutine X0_find(this, metr, x, step, x0)

		Class(interp) :: this
		Class(metric) :: metr
		Integer(4), intent(in) :: x, step
		Integer(4), intent(out) :: x0
		Integer(4) :: xk
		Real(8) :: gap

		do xk = 1, 2*this.dim

			gap = (metr.latlon_c(1, 1 - step, x, 2)) - (metr.latlon_c(1, step, xk, 2))
			if(gap >= 0d0) then
				x0 = xk
			end if

		end do

	end Subroutine



	Subroutine weight_find(this, metr, x0, x, step, weight)

		Class(interp) :: this
		Class(metric) :: metr
		Type(geometry) :: geom
		Integer(4), intent(in) :: step, x, x0
		Real(8), intent(out) :: weight(1:this.n)
		Integer(4) :: n, y(1:this.n), i, j, dim
		Real(8) :: s, numen, k(4, 4), h(1:4), latlon_y(this.n), latlon_x

		n = this.n;  dim = this.dim; weight(:) = 1d0
		latlon_x = metr.latlon_c(1,1-step,x,2)

		y(1) = x0 - n/2 + 1
		latlon_y(1) = metr.latlon_c(1,2*dim+1-step,y(1),5)
		h(1) = abs(latlon_x - latlon_y(1))

		do i = 2, n
			y(i) = y(i-1)+1
			latlon_y(i) = metr.latlon_c(1,2*dim+1-step,y(i),5)
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


		latlon_x = metr.latlon_c(1,2*dim+1,1,5)

		do i = 1, 2
			y(i) = i
			latlon_y(i) = metr.latlon_c(1,1,y(i),2)
			h(i) = latlon_x - latlon_y(i)
			do j = 1, 2
				if(j /= i) k(i, j) = h(i)/(latlon_y(j) - latlon_y(i))
				if(j /= i) k(j, i) = h(j)/(latlon_y(i) - latlon_y(j))
			end do
			this.fantom1(i) = 1d0
		end do

		do i = 1, 2
			do j = 1, 2
				if(j /= i) this.fantom1(i) = this.fantom1(i)*k(j,i)
			end do
		end do

	end Subroutine



end module