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
		Class(f_var) :: var, var_pr
		Class(g_var) :: grid

		real(8) g, height, dt, partial(1:2)
		integer(4) face_idx, x, y, dim, i, j

		g = grid.g;  height = var_pr.height;  dim = var_pr.dim

		dt = 10000d0


		do face_idx = 1, 6
			do y = -dim, dim
				do x = -dim, dim

					partial(1) = grid.partial_c1_x(var_pr.h_height(x-1:x+1, y, face_idx), x, y)
					var.u_vel(x, y, face_idx) = var_pr.u_vel(x, y, face_idx) - dt*g*partial(1)

					partial(1) = grid.partial_c1_y(var_pr.h_height(x, y-1:y+1, face_idx), x, y)
					var.v_vel(x, y, face_idx) = var_pr.v_vel(x, y, face_idx) - dt*g*partial(1)

					partial(1) = grid.partial_c1_x(var_pr.u_vel(x-1:x+1, y, face_idx), x, y)
					partial(2) = grid.partial_c1_y(var_pr.v_vel(x, y-1:y+1, face_idx), x, y)
					var.h_height(x, y, face_idx) = var_pr.h_height(x, y, face_idx) - height*(partial(1) + partial(2))


		! print '(" step = ", f10.2)', partial(1)
		! print '(" step = ", f10.2)', partial(1)


				end do
			end do

			do y = -dim, dim, 2*dim+1
				do x = -dim, dim, 2*dim+1
					do i = 1, var.step
						do j = 1, var.step

					var.u_vel(x, y, face_idx) = (var_pr.u_vel(x+i, y+j, face_idx) + var_pr.u_vel(x-i, y+j, face_idx) + var_pr.u_vel(x+i, y-j, face_idx) + var_pr.u_vel(x-i, y-j, face_idx))/4
					var.v_vel(x, y, face_idx) = (var_pr.v_vel(x+i, y+j, face_idx) + var_pr.v_vel(x-i, y+j, face_idx) + var_pr.v_vel(x+i, y-j, face_idx) + var_pr.v_vel(x-i, y-j, face_idx))/4
					var.h_height(x, y, face_idx) = (var_pr.h_height(x+i, y+j, face_idx) + var_pr.h_height(x-i, y+j, face_idx) + var_pr.h_height(x+i, y-j, face_idx) + var_pr.h_height(x-i, y-j, face_idx))/4

						end do
					end do
				end do
			end do

		end do

		call var_pr.equal(var)
		call var_pr.borders(var_pr)
		call var_pr.corner_zero()


	end subroutine




end module
