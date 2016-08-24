module derivatives

	use grid_var, Only: g_var

implicit none

	Private
	Public :: der

	Type der

		CONTAINS
		Procedure, Public :: div_2 => div_2
		Procedure, Public :: div_4 => div_4
		Procedure, Public :: partial_c2_x => partial_c2_x
		Procedure, Public :: partial_c2_y => partial_c2_y
		Procedure, Public :: partial_c4_x => partial_c4_x
		Procedure, Public :: partial_c4_y => partial_c4_y
	End Type


CONTAINS


	real(8) function partial_c2_x(this, g, fun, x, y)
		Class(der) :: this
		Class(g_var) :: g
		real(8), intent(in) :: fun(-1:1)
		integer, intent(in) :: x, y
		partial_c2_x = (fun(1) * g.x_dist(x, y) - fun(-1) * g.x_dist(x+1, y))/&
		(g.x_dist(x, y) + g.x_dist(x+1, y))**2

	end function



	real(8) function partial_c2_y(this, g, fun, x, y)
		Class(der) :: this
		Class(g_var) :: g
		real(8), intent(in) :: fun(-1:1)
		integer, intent(in) :: x, y
		partial_c2_y = (fun(1) * g.y_dist(x, y) - fun(-1)*g.y_dist(x, y+1))/&
		(g.y_dist(x, y) + g.y_dist(x, y+1))**2

	end function



	real(8) function div_2(this, g, u1, u2, x, y)
		Class(der) :: this
		Class(g_var) :: g
		real(8), intent(in) :: u1(-1:1), u2(-1:1)
		integer(4), intent(in) :: x, y
		integer(4) i
		real(8) u_1(-1:1), u_2(-1:1), G_11, G_12, G_21, G_22

		do i = -1, 1
			G_11 = g.G_inverse(x+i, y, 1, 1)
			G_12 = g.G_inverse(x, y+i, 1, 2)
			G_21 = g.G_inverse(x+i, y, 2, 1)
			G_22 = g.G_inverse(x, y+i, 2, 2)
			u_1(i) = G_11*u1(i)! + G_12*u2(i)
			u_2(i) = G_22*u2(i)! + G_21*u1(i)
		end do



		div_2 = ((u_1(1) * g.x_dist(x, y)*g.G_sqr(x+1,y) - u_1(-1)*g.x_dist(x+1, y)*g.G_sqr(x-1,y))/&
			(g.x_dist(x, y) + g.x_dist(x+1, y))**2 &
		+ (u_2(1) * g.y_dist(x, y)*g.G_sqr(x,y+1) - u_2(-1)*g.y_dist(x, y+1)*g.G_sqr(x,y-1))/&
			(g.y_dist(x, y) + g.y_dist(x, y+1))**2)/g.G_sqr(x,y)

	end function



	real(8) function div_4(this, g, u1, u2, x, y)
		Class(der) :: this
		Class(g_var) :: g
		real(8), intent(in) :: u1(-2:2), u2(-2:2)
		integer(4), intent(in) :: x, y
		integer(4) i
		real(8) u_1(-2:2), u_2(-2:2), G_11, G_12, G_21, G_22, J_1(-2:2), J_2(-2:2)
		real(8) Ax , Bx, Cx, Dx, Ex, Ay, By, Cy, Dy, Ey

		Ax = g.four_order_const_x( 1, x, y);  Bx = g.four_order_const_x( 2, x, y);  Cx = g.four_order_const_x( 3, x, y)
		Dx = g.four_order_const_x( 4, x, y);  Ex = g.four_order_const_x( 5, x, y)

		Ay = g.four_order_const_y( 1, x, y);  By = g.four_order_const_y( 2, x, y);  Cy = g.four_order_const_y( 3, x, y)
		Dy = g.four_order_const_y( 4, x, y);  Ey = g.four_order_const_y( 5, x, y)

		do i = -2, 2
			G_11 = g.G_inverse(x+i, y, 1, 1)
			G_12 = g.G_inverse(x, y+i, 1, 2)
			G_21 = g.G_inverse(x+i, y, 2, 1)
			G_22 = g.G_inverse(x, y+i, 2, 2)
			u_1(i) = G_11*u1(i)! + G_12*u2(i)
			u_2(i) = G_22*u2(i)! + G_21*u1(i)
			J_1(i) = g.G_sqr(x+i, y)
			J_2(i) = g.G_sqr(x, y+i)
		end do


		div_4 = (Ax*u_1(1)*J_1(1) + Bx*u_1(-1)*J_1(-1) + Cx*u_1(2)*J_1(2) + Dx*u_1(-2)*J_1(-2) + Ex*u_1(0)*J_1(0) + &
			 Ay*u_2(1)*J_2(1) + By*u_2(-1)*J_2(-1) + Cy*u_2(2)*J_2(2) + Dy*u_2(-2)*J_2(-2) + Ey*u_2(0)*J_2(0))/J_1(0)

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