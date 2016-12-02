module derivatives

	use metrics, Only: metric

implicit none

	Private
	Public :: der

	Type der

		CONTAINS
		Procedure, Public :: partial_c => partial_c
		Procedure, Private :: partial_c2 => partial_c2
		Procedure, Private :: partial_c4 => partial_c4
		Procedure, Private :: partial_c6 => partial_c6
		Procedure, Private :: partial_c_fg => partial_c_fg
		Procedure, Public :: div => div
	End Type


CONTAINS


	real(8) function partial_c2(this, fun, h)
		Class(der) :: this
		real(8), intent(in) :: fun(-1:1), h
		partial_c2 = (fun(1) - fun(-1))/(2d0*h)

	end function



	real(8) function partial_c4(this, fun, h)
		Class(der) :: this
		real(8), intent(in) :: fun(-2:2), h
		real(8) A , B, C, D

		A = 2d0/(3d0*h);  B = - 2d0/(3d0*h);  C = - 1.0/(12d0*h);  D = 1d0/(12d0*h)
		partial_c4 = A*fun(1) + B*fun(-1) + C*fun(2) + D*fun(-2)

	end function


	real(8) function partial_c6(this, fun, h)
		Class(der) :: this
		real(8), intent(in) :: fun(-3:3), h
		real(8) A , B, C, D, E, F

		A = 3d0/(4d0*h);  B = - 3d0/(4d0*h);  C = - 3d0/(20d0*h);  D = 3d0/(20d0*h);  E = 1d0/(60d0*h);  F = - 1d0/(60d0*h)
		partial_c6 = A*fun(1) + B*fun(-1) + C*fun(2) + D*fun(-2) + E*fun(3) + F*fun(-3)

	end function



	real(8) function partial_c(this, fun, h, step)
		Class(der) :: this
		integer(4), intent(in) :: step
		real(8), intent(in) :: fun(-step:step), h

		if (step == 1) then
			partial_c = this.partial_c2(fun, h)
		else if (step == 2) then
			partial_c = this.partial_c4(fun, h)
		else if (step == 3) then
			partial_c = this.partial_c6(fun, h)
		end if

	end function



	real(8) function partial_c_fg(this, f_fun, g_fun, h, step)
		Class(der) :: this
		integer(4), intent(in) :: step
		real(8), intent(in) :: f_fun(-step:step), g_fun(-step:step), h

			partial_c_fg = g_fun(0)*this.partial_c(f_fun, h, step) + f_fun(0)*this.partial_c(g_fun, h, step)

	end function



	real(8) function div(this, metr, u1_con, u2_con, h, x, y, step)
		Class(der) :: this
		Class(metric) :: metr
		integer(4), intent(in) :: x, y, step
		real(8), intent(in) :: u1_con(-step:step), u2_con(-step:step), h
		integer(4) i, j
		real(8) J_1(-step:step), J_2(-step:step), G(2,2)

		do i = -step, step
			J_1(i) = metr.G_sqr(x+i, y)
			J_2(i) = metr.G_sqr(x, y+i)
		end do

		div = ( this.partial_c_fg(u1_con, J_1, h, step) + this.partial_c_fg(u2_con, J_2, h, step))/J_1(0)

	end function


end module