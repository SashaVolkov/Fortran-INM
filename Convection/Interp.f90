Module Interp

	Use modnet

	IMPLICIT NONE

	Public :: intp

	Type intp


	CONTAINS

! 		Procedure :: init=> intp_init
		Procedure :: Lintp=> intp_Lintp

	End Type

	CONTAINS

	Subroutine intp_Lintp(this, x_wave, Mass_in, res, g)

		Class(grid) :: g
		Class(intp) :: this

		Real(8), Intent(in) :: Mass_in(g.ns - g.bstep : g.nf + g.fstep)
		Real(8), Intent(in) :: x_wave
		Real(8), Intent(out) :: res
		Integer(4) xj, xk, k, j, n, x0, set_x
		Real(8) :: temp1, temp2, res_min


		n=1
		set_x = floor(x_wave)
		x0 = set_x - floor(((n*1.0)/2.0))
		xk = x0;  temp2=0.0;  temp1=1.0

		do k=1,n+1
			xj=x0

			do j = 1,n+1
				if ( xj /= xk ) temp1 = temp1*((x_wave - xj)*1.0/(xk-xj)*1.0)
				xj=xj+1
			end do

			temp2 = temp2 + (Mass_in(xk)*temp1)
			temp1=1.0
			xk=xk+1

		end do

		res = temp2

	End Subroutine

End Module