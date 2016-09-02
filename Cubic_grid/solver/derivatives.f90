module derivatives

	use metrics, Only: metric
	use grid_var, Only: g_var

implicit none

	Private
	Public :: der

	Type der

		CONTAINS
		Procedure, Public :: div_2 => div_2
		Procedure, Public :: div_4 => div_4
		Procedure, Public :: partial_c2 => partial_c2
		Procedure, Public :: partial_c4_x => partial_c4_x
		Procedure, Public :: partial_c4_y => partial_c4_y
	End Type


CONTAINS


	real(8) function partial_c2(this, fun, dist)
		Class(der) :: this
		real(8), intent(in) :: fun(-1:1), dist(0:1)
		partial_c2 = (fun(1) - fun(-1))/(dist(0) + dist(1))

	end function


	real(8) function div_2(this, grid, metr, u1_cov, u2_cov, x, y)
		Class(der) :: this
		Class(g_var) :: grid
		Class(metric) :: metr
		real(8), intent(in) :: u1_cov(-1:1), u2_cov(-1:1)
		integer(4), intent(in) :: x, y
		integer(4) i, j
		real(8) u1_con(-1:1), u2_con(-1:1), J_1(-1:1), J_2(-1:1), G(2,2), dx, dy

		do i = -1, 1
			J_1(i) = metr.G_sqr(x+i, y)
			J_2(i) = metr.G_sqr(x, y+i)
			G(1,1) = metr.G_inverse(1, 1, x+i, y)
			G(1,2) = metr.G_inverse(1, 2, x, y+i)
			G(2,1) = metr.G_inverse(2, 1, x+i, y)
			G(2,2) = metr.G_inverse(2, 2, x, y+i)
			u1_con(i) = G(1,1)*u1_cov(i) + G(1,2)*u2_cov(i)
			u2_con(i) = G(2,2)*u2_cov(i) + G(2,1)*u1_cov(i)
			dx = 2*grid.delta_on_cube !grid.x_dist(x, y) + grid.x_dist(x+1, y)
			dy = 2*grid.delta_on_cube !grid.y_dist(x, y) + grid.y_dist(x, y+1)
		end do


		div_2 = ( (u1_con(1)*J_1(1) - u1_con(-1)*J_1(-1))/dx + (u2_con(1)*J_2(1) - u2_con(-1)*J_2(-1))/dy )/J_1(0)

	end function




	real(8) function div_4(this, g, metr, u1_cov, u2_cov, x, y)
		Class(der) :: this
		Class(g_var) :: g
		Class(metric) :: metr
		real(8), intent(in) :: u1_cov(-2:2), u2_cov(-2:2)
		integer(4), intent(in) :: x, y
		integer(4) i
		real(8) u1_con(-2:2), u2_con(-2:2), G_11, G_12, G_21, G_22, J_1(-2:2), J_2(-2:2)
		real(8) Ax , Bx, Cx, Dx, Ex, Ay, By, Cy, Dy, Ey

		Ax = g.four_order_const_x( 1, x, y);  Bx = g.four_order_const_x( 2, x, y);  Cx = g.four_order_const_x( 3, x, y)
		Dx = g.four_order_const_x( 4, x, y);  Ex = g.four_order_const_x( 5, x, y)

		Ay = g.four_order_const_y( 1, x, y);  By = g.four_order_const_y( 2, x, y);  Cy = g.four_order_const_y( 3, x, y)
		Dy = g.four_order_const_y( 4, x, y);  Ey = g.four_order_const_y( 5, x, y)

		do i = -2, 2
			G_11 = metr.G_inverse(1, 1, x+i, y)
			G_12 = metr.G_inverse(1, 2, x, y+i)
			G_21 = metr.G_inverse(2, 1, x+i, y)
			G_22 = metr.G_inverse(2, 2, x, y+i)
			u1_con(i) = G_11*u1_cov(i) + G_12*u2_cov(i)
			u2_con(i) = G_22*u2_cov(i) + G_21*u1_cov(i)
			J_1(i) = metr.G_sqr(x+i, y)
			J_2(i) = metr.G_sqr(x, y+i)
		end do


		div_4 = ( Ax*u1_con(1)*J_1(1) + Bx*u1_con(-1)*J_1(-1) + Cx*u1_con(2)*J_1(2) + Dx*u1_con(-2)*J_1(-2) + Ex*u1_con(0)*J_1(0) + &
			 Ay*u2_con(1)*J_2(1) + By*u2_con(-1)*J_2(-1) + Cy*u2_con(2)*J_2(2) + Dy*u2_con(-2)*J_2(-2) + Ey*u2_con(0)*J_2(0) )/J_1(0)

	end function



	real(8) function partial_c4_x(this, g, fun, x, y)
		Class(der) :: this
		Class(g_var) :: g
		real(8), intent(in) :: fun(-2:2)
		integer, intent(in) :: x, y
		real(8) A , B, C, D, E

		A = g.four_order_const_x( 1, x, y);  B = g.four_order_const_x( 2, x, y);  C = g.four_order_const_x( 3, x, y)
		D = g.four_order_const_x( 4, x, y);  E = g.four_order_const_x( 5, x, y)

		partial_c4_x = A*fun(1) + B*fun(-1) + C*fun(2) + D*fun(-2) +  E*fun(0)

	end function


	real(8) function partial_c4_y(this, g, fun, x, y)
		Class(der) :: this
		Class(g_var) :: g
		real(8), intent(in) :: fun(-2:2)
		integer, intent(in) :: x, y
		real(8) A , B, C, D, E

		A = g.four_order_const_y( 1, x, y);  B = g.four_order_const_y( 2, x, y);  C = g.four_order_const_y( 3, x, y)
		D = g.four_order_const_y( 4, x, y);  E = g.four_order_const_y( 5, x, y)

		partial_c4_y = A*fun(1) + B*fun(-1) + C*fun(2) + D*fun(-2) +  E*fun(0)

	end function



end module