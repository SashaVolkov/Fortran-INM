module interpolation

	use grid_var, Only: g_var
	use sphere_geometry, Only: geometry
	use mpi

	implicit none


	Private
	Public :: interp

	Type interp

		Real(8), Allocatable :: weight(:, :, :)
		Integer(8), Allocatable :: x0_mass(:, :)
		integer(4) ns, nf, n

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
		Integer(4) :: k, i, j, xk, x, x0

		this.ns = g.ns_xy(1);  this.nf = g.nf_xy(1);  this.n = n

		Allocate(this.x0_mass(this.ns: this.nf, 1: g.step))
		Allocate(this.weight(this.ns: this.nf, this.ns: this.nf, 1: g.step))

		do i = 1, g.step
			do x = this.ns, this.nf
				call this.X0_find(g, x, i, x0)
				this.x0_mass(x, i) = x0
				print *,  this.x0_mass(x, i), x
			end do
		end do


		do i = 1, g.step
			do x = this.ns, this.nf
				do xk = this.x0_mass(x, i) - n/2, this.x0_mass(x, i) + (n+1)/2
					this.weight(xk, x, i) = this.weight_find(g, xk, x, i)
				end do
			end do
		end do

	End Subroutine


	Subroutine Lagrange(this, Mass, step)

		Class(interp) :: this

		Real(8), Intent(inout) :: Mass(this.ns: this.nf)
		Integer(4), intent(in) :: step
		Integer(4) xj, xk, k, j, n, x0, set_x, x, s
		Real(8) :: temp, Mass_temp(this.ns: this.nf)

		n=this.n

		do x = this.ns, this.nf
			x0 = this.x0_mass(x, step)
			s = x0 - n/2
			temp=0.0

			do xk=s, s + n
				temp = temp + (Mass(xk)*this.weight(xk, x, step))
			end do

			Mass_temp(x) = temp
		end do

		Mass(:) = Mass_temp(:)

	End Subroutine



	Subroutine X0_find(this, g, x, step, x0)

		Class(interp) :: this
		Class(g_var) :: g
		Integer(4), intent(in) :: x, step
		Integer(4), intent(out) :: x0
		Integer(4) :: xk
		Real(8) :: gap

		do xk = this.ns, this.nf

			gap = (g.latlon_c(1, 1 - step, x, 2)) - (g.latlon_c(1, step, xk, 2))
			! print *, g.latlon_c(1, 1 - step, x, 2), g.latlon_c(1, step, xk, 2), gap
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
					temp = (geom.dist(latlon_x(:), latlon_xk(:)))/(geom.dist(latlon_xj(:), latlon_xk(:)))
				end if
			end do

		weight_find = temp


	end function



end module
