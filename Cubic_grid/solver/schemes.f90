module schemes

	use grid_var, Only: g_var
	use func_var, Only: f_var

	implicit none

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
		integer(4) face, x, y, dim, i, j, stat

		g = grid.g;  height = var_pr(1).height;  dim = var_pr(1).dim
		dt = grid.dt


		do face = 1, 6
			do y = -dim, dim
				do x = -dim, dim

					partial(1) = grid.partial_c1_x(var_pr(face).h_height(x-1:x+1, y), x, y)
					var(face).u_vel(x, y) = var_pr(face).u_vel(x, y) - dt*g*partial(1)

					partial(1) = grid.partial_c1_y(var_pr(face).h_height(x, y-1:y+1), x, y)
					var(face).v_vel(x, y) = var_pr(face).v_vel(x, y) - dt*g*partial(1)

					partial(1) = grid.partial_c1_x(var_pr(face).u_vel(x-1:x+1, y), x, y)
					partial(2) = grid.partial_c1_y(var_pr(face).v_vel(x, y-1:y+1), x, y)
					var(face).h_height(x, y) = var_pr(face).h_height(x, y) - height*(partial(1) + partial(2))


				end do
			end do

			do y = -dim, dim, 2*dim+1
				do x = -dim, dim, 2*dim+1
					do i = 1, var(1).step
						do j = 1, var(1).step

					var(face).u_vel(x, y) = (var_pr(face).u_vel(x+i, y+j) + var_pr(face).u_vel(x-i, y+j) +&
					 var_pr(face).u_vel(x+i, y-j) + var_pr(face).u_vel(x-i, y-j))/4

					var(face).v_vel(x, y) = (var_pr(face).v_vel(x+i, y+j) + var_pr(face).v_vel(x-i, y+j) +&
					 var_pr(face).v_vel(x+i, y-j) + var_pr(face).v_vel(x-i, y-j))/4

					var(face).h_height(x, y) = (var_pr(face).h_height(x+i, y+j) + var_pr(face).h_height(x-i, y+j) +&
					 var_pr(face).h_height(x+i, y-j) + var_pr(face).h_height(x-i, y-j))/4

						end do
					end do
				end do
			end do


		call var_pr(face).equal(var(face))

		end do
		! stat = var(face).borders(var, var_pr)



	end subroutine




end module
