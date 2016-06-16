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
		Class(f_var) :: var(1:6), var_pr(1:6)
		Class(g_var) :: grid

		real(8) g, height, dt, partial(1:2)
		integer(4) face, x, y, dim, i, j, stat, ns_x, ns_y, nf_x, nf_y, ier

		g = grid.g;  height = var_pr(1).height;  dim = var_pr(1).dim
		dt = grid.dt;  ns_x = var(1).ns_x;  ns_y = var(1).ns_y
		nf_x = var(1).nf_x;  nf_y = var(1).nf_y


		do face = 1, 6
			do y = ns_y, nf_y
				do x = ns_x, nf_x

					partial(1) = grid.partial_c1_x(var_pr(face).h_height(x-1:x+1, y), x, y)
					var(face).u_vel(x, y) = var_pr(face).u_vel(x, y) - dt*g*partial(1)

					partial(1) = grid.partial_c1_y(var_pr(face).h_height(x, y-1:y+1), x, y)
					var(face).v_vel(x, y) = var_pr(face).v_vel(x, y) - dt*g*partial(1)

					partial(1) = grid.partial_c1_x(var_pr(face).u_vel(x-1:x+1, y), x, y)
					partial(2) = grid.partial_c1_y(var_pr(face).v_vel(x, y-1:y+1), x, y)
					var(face).h_height(x, y) = var_pr(face).h_height(x, y) - height*(partial(1) + partial(2))


				end do
			end do
		end do

		do face = 1, 6
			call var_pr(face).equal(var(face))
		end do


	end subroutine




end module
