module derivatives

	use metrics, Only: metric

implicit none

	Private
	Public :: der

	Type der

		CONTAINS
		Procedure, Public :: partial_c => partial_c
		Procedure, Public :: partial_c2 => partial_c2
		Procedure, Public :: partial_c_fg => partial_c_fg
		Procedure, Public :: div => div
		Procedure, Public :: curl => curl
		Procedure, Public :: grad_div_vec => grad_div_vec
		Procedure, Public :: laplace => laplace
	End Type


CONTAINS



	Real(8) function partial_c(this, fun, h, step)
		Class(der) :: this
		Integer(4), intent(in) :: step
		Real(8), intent(in) :: fun(-step:step), h
		Real(8) A , B, C, D, E, F, G, I

		if (step == 1) then
			partial_c = (fun(1) - fun(-1))/(2d0*h)
		else if(step == 2) then
			A = 2d0/(3d0*h);  B = - 2d0/(3d0*h);  C = - 1d0/(12d0*h);  D = 1d0/(12d0*h)
			partial_c = A*fun(1) + B*fun(-1) + C*fun(2) + D*fun(-2)
		else if(step == 3) then
			A = 3d0/(4d0*h);  B = - 3d0/(4d0*h);  C = - 3d0/(20d0*h);  D = 3d0/(20d0*h)
			E = 1d0/(60d0*h);  F = - 1d0/(60d0*h)
			partial_c = A*fun(1) + B*fun(-1) + C*fun(2) + D*fun(-2) + E*fun(3) + F*fun(-3)
		else if(step == 4) then
			A = 4d0/(5d0*h);  B = - 4d0/(5d0*h);  C = - 1d0/(5d0*h);  D = 1d0/(5d0*h)
			E = 4d0/(105d0*h);  F = - 4d0/(105d0*h); G = - 1d0/(280d0*h); I = 1d0/(280d0*h)

			partial_c = A*fun(1) + B*fun(-1) + C*fun(2) + D*fun(-2) + E*fun(3) + F*fun(-3) + G*fun(4) + I*fun(-4)
		end if

	end function



	Real(8) function partial_c2(this, fun, h, step)
		Class(der) :: this
		Integer(4), intent(in) :: step
		Real(8), intent(in) :: fun(-step:step), h
		Real(8) A , B, C, D, E, F, G, I, J

		if (step == 1) then
			partial_c2 = (fun(1) + fun(-1) - 2d0*fun(0))/(h)
		else if(step == 2) then
			A = 4d0/(3d0*h);  B = 4d0/(3d0*h);  C = - 1d0/(12d0*h);  D = -1d0/(12d0*h);  E = -5d0/(2d0*h)
			partial_c2 = A*fun(1) + B*fun(-1) + C*fun(2) + D*fun(-2) + E*fun(0)
		else if(step == 3) then
			A = 3d0/(2d0*h);  B = 3d0/(2d0*h);  C = - 3d0/(20d0*h);  D = -3d0/(20d0*h)
			E = 1d0/(90d0*h);  F = 1d0/(90d0*h); G = -49d0/(18d0*h)
			partial_c2 = A*fun(1) + B*fun(-1) + C*fun(2) + D*fun(-2) + E*fun(3) + F*fun(-3) + G*fun(0)
		else if(step == 4) then
			A = 8d0/(5d0*h);  B = 8d0/(5d0*h);  C = - 1d0/(5d0*h);  D = - 1d0/(5d0*h)
			E = 8d0/(315d0*h);  F = 8d0/(315d0*h); G = - 1d0/(560d0*h); I = -1d0/(560d0*h);  J = -205d0/(72d0)

			partial_c2 = A*fun(1) + B*fun(-1) + C*fun(2) + D*fun(-2) + E*fun(3) + F*fun(-3) + G*fun(4) + I*fun(-4) + J*fun(0)
		end if

	end function



	Real(8) function partial_c_fg(this, f_fun, g_fun, h, step)
		Class(der) :: this
		Integer(4), intent(in) :: step
		Real(8), intent(in) :: f_fun(-step:step), g_fun(-step:step), h

			partial_c_fg = g_fun(0)*this.partial_c(f_fun, h, step) + f_fun(0)*this.partial_c(g_fun, h, step)

	end function



	Real(8) function div(this, metr, u1_con, u2_con, h, x, y, step)
		Class(der) :: this
		Class(metric) :: metr
		Integer(4), intent(in) :: x, y, step
		Real(8), intent(in) :: u1_con(-step:step), u2_con(-step:step), h
		Integer(4) i, j
		Real(8) J_1(-step:step), J_2(-step:step), G(2,2)

		do i = -step, step
			J_1(i) = metr.G_sqr(x+i, y)
			J_2(i) = metr.G_sqr(x, y+i)
		end do

		div = ( this.partial_c_fg(u1_con, J_1, h, step) + this.partial_c_fg(u2_con, J_2, h, step))/J_1(0)

	end function



	Real(8) function curl(this, metr, u1_con, u2_con, h, x, y, step)
		Class(der) :: this
		Class(metric) :: metr
		Integer(4), intent(in) :: x, y, step
		Real(8), intent(in) :: u1_con(-step:step), u2_con(-step:step), h
		Real(8) :: u1_cov(-step:step), u2_cov(-step:step)
		Integer :: i

		do i = -step, step
			call metr.con_to_cov(u1_con(i),u2_con(i),u1_cov(i),u2_cov(i),x,y)
		end do

		curl = ( this.partial_c(u2_cov, h, step) - this.partial_c(u1_cov, h, step))/metr.G_sqr(x, y)

	end function


	Subroutine grad_div_vec(this, metr, u_con, v_con, u, v, h, x, y, step)
		Class(der) :: this
		Class(metric) :: metr
		Integer(4), intent(in) :: x, y, step
		Real(8), intent(in) :: u_con(-step:step, -step:step), v_con(-step:step, -step:step), h
		Real(8), intent(out) :: u, v
		Real(8) J_1(-step:step), J_2(-step:step), div, temp_u(-step:step), temp_v(-step:step)
		Integer :: i

		temp_v = 1d0;  temp_u = 1d0

		do i = -step, step
			J_1(i) = metr.G_sqr(x+i, y)
			J_2(i) = metr.G_sqr(x, y+i)
			temp_u(i) = this.partial_c(metr.G_sqr(x-step:x+step, y+i)*v_con(:,i), h, step)
			temp_v(i) = this.partial_c(metr.G_sqr(x-step:x+step, y+i)*u_con(:,i), h, step)
		end do
		div = this.div(metr, u_con, v_con, h, x, y, step)

		u = div*this.partial_c(1d0/J_1, h, step) + this.partial_c2(J_1*u_con(:, 0), h, step)/J_1(0) + this.partial_c(temp_u, h, step)/J_1(0)
		v = div*this.partial_c(1d0/J_2, h, step) + this.partial_c2(J_1*v_con(0, :), h, step)/J_1(0) + this.partial_c(temp_v, h, step)/J_1(0)

	end Subroutine


	Real(8) function laplace(this, metr, f, h, x, y, step)
		Class(der) :: this
		Class(metric) :: metr
		Integer(4), intent(in) :: x, y, step
		Real(8), intent(in) :: f(-step:step, -step:step), h
		Real(8) :: f_x(-step:step), f_y(-step:step), sum_x(-step:step), sum_y(-step:step)
		Real(8) J_1(-step:step), J_2(-step:step)
		Integer :: i

		do i = -step, step
			f_x = f(-step:step, i)
			f_y = f(i, -step:step)

			J_1(i) = metr.G_sqr(x+i, y)
			J_2(i) = metr.G_sqr(x, y+i)

			sum_x(i) = J_1(i)*( metr.G_tensor(1, 1, x, y+i)*this.partial_c(f_x, h, step) + metr.G_tensor(1, 2, x+i, y)*this.partial_c(f_y, h, step))
			sum_y(i) = J_2(i)*( metr.G_tensor(2, 1, x, y+i)*this.partial_c(f_x, h, step) + metr.G_tensor(2, 2, x+i, y)*this.partial_c(f_y, h, step))
		end do

		laplace = ( this.partial_c(sum_x, h, step) + this.partial_c(sum_y, h, step))/J_1(0)

	end function


end module