module schemes

	use grid_var, Only: g_var
	use func_var, Only: f_var

	implicit none

	include"mpif.h"

	Private
	Public :: schema

	Type schema
		CONTAINS
		Procedure, Public :: Linear => Linear
	End Type

	CONTAINS



	subroutine Linear(this, var, var_pr, grid)

		Class(schema) :: this
		Class(f_var) :: var, var_pr
		Class(g_var) :: grid

		real(8) g, height, dt, partial(1:2)
		integer(4) face, x, y, dim, i, j, stat, ns_x, ns_y, nf_x, nf_y, ier

		g = grid.g;  height = var_pr.height;  dim = var_pr.dim
		dt = grid.dt;  ns_x = var.ns_x;  ns_y = var.ns_y
		nf_x = var.nf_x;  nf_y = var.nf_y


		do face = 1, 6
			do y = ns_y, nf_y
				do x = ns_x, nf_x

					partial(1) = grid.partial_c1_x(var_pr.h_height(x-1:x+1, y, face), x, y)
					var.u_vel(x, y, face) = var_pr.u_vel(x, y, face) - dt*g*partial(1)

					partial(1) = grid.partial_c1_y(var_pr.h_height(x, y-1:y+1, face), x, y)
					var.v_vel(x, y, face) = var_pr.v_vel(x, y, face) - dt*g*partial(1)

					partial(1) = grid.partial_c1_x(var_pr.u_vel(x-1:x+1, y, face), x, y)
					partial(2) = grid.partial_c1_y(var_pr.v_vel(x, y-1:y+1, face), x, y)
					var.h_height(x, y, face) = var_pr.h_height(x, y, face) - height*(partial(1) + partial(2))


				end do
			end do
		end do

			call var_pr.equal(var)



	end subroutine




end module
