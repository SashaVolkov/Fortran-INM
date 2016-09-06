module derivatives

	use metrics, Only: metric

implicit none

	Private
	Public :: der

	Type der

		CONTAINS
		Procedure, Public :: div => div
		Procedure, Public :: partial_c2 => partial_c2
		! Procedure, Public :: partial_c2_fg => partial_c2_fg
		Procedure, Public :: partial_c4 => partial_c4
		Procedure, Public :: partial_c_fg => partial_c_fg
	End Type


CONTAINS


	real(8) function partial_c2(this, fun, h)
		Class(der) :: this
		real(8), intent(in) :: fun(-1:1), h
		partial_c2 = (fun(1) - fun(-1))/(2.0*h)

	end function



	real(8) function partial_c4(this, fun, h)
		Class(der) :: this
		real(8), intent(in) :: fun(-2:2), h
		real(8) A , B, C, D, E

		A = 2.0/(3.0*h);  B = - 2.0/(3.0*h);  C = - 1.0/(12.0*h);  D = 1.0/(12.0*h);  E = 0.0
		partial_c4 = A*fun(1) + B*fun(-1) + C*fun(2) + D*fun(-2) +  E*fun(0)

	end function



	real(8) function partial_c_fg(this, f_fun, g_fun, h, order)
		Class(der) :: this
		integer(4), intent(in) :: order
		real(8), intent(in) :: f_fun(-order:order), g_fun(-order:order), h

		if (order == 1) then
			partial_c_fg = g_fun(0)*this.partial_c2(f_fun, h) + f_fun(0)*this.partial_c2(g_fun, h)
		else if (order == 2) then
			partial_c_fg = g_fun(0)*this.partial_c4(f_fun, h) + f_fun(0)*this.partial_c4(g_fun, h)
		end if

	end function



	real(8) function div(this, metr, u1_cov, u2_cov, h, x, y, order)
		Class(der) :: this
		Class(metric) :: metr
		integer(4), intent(in) :: x, y, order
		real(8), intent(in) :: u1_cov(-order:order), u2_cov(-order:order), h
		integer(4) i
		real(8) u1_con(-order:order), u2_con(-order:order), J_1(-order:order), J_2(-order:order), G(2,2)

		do i = -order, order
			J_1(i) = metr.G_sqr(x+i, y)
			J_2(i) = metr.G_sqr(x, y+i)
			G(1,1) = metr.G_inverse(1, 1, x+i, y)
			G(1,2) = metr.G_inverse(1, 2, x+i, y)
			G(2,1) = metr.G_inverse(2, 1, x, y+i)
			G(2,2) = metr.G_inverse(2, 2, x, y+i)
			u1_con(i) = G(1,1)*u1_cov(i) + G(1,2)*u2_cov(i)
			u2_con(i) = G(2,2)*u2_cov(i) + G(2,1)*u1_cov(i)
		end do

		div = ( this.partial_c_fg(u1_con, J_1, h, order) + this.partial_c_fg(u2_con, J_2, h, order))/J_1(0)

	end function


end module