module interpolation

	use grid_var, Only: g_var
	use func_var, Only: f_var
	use mpi

	implicit none


	Private
	Public :: interp

	Type interp

		Real(8), Allocatable :: weight(:, :, :)
		integer(4) ns, nf

		CONTAINS
! 		Procedure, Public :: init => init
		Procedure, Public :: Lagrange => Lagrange
	End Type

	CONTAINS



	Subroutine Lagrange(this, Mass)

		Class(interp) :: this

		Real(8), Intent(inout) :: Mass(this.ns: this.nf)
		Integer(4) xj, xk, k, j, n, x0, set_x, x
		Real(8) :: temp, Mass_temp(this.ns: this.nf)

		n=1

		do x = this.ns, this.nf
			x0 = x - floor(((n*1.0)/2.0))
			temp=0.0

			do xk=x0, x0 + n
				temp = temp + (Mass(xk)*this.weight(xk, x0, 1))
			end do

			Mass_temp(x) = temp
		end do

		Mass(:) = Mass_temp(:)

	End Subroutine

end modul1e
