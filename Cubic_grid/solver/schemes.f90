module schemes

	use special_variables, Only: variables
	use var_func, Only: func

	implicit none

	Private
	Public :: schema

	Type schema
		CONTAINS
		Procedure, Public :: Linear => Linear
	End Type

	CONTAINS



	subroutine Linear(this, var, var_pr, f)

		Class(schema) :: this
		Class(variables) :: var, var_pr
		Class(func) :: f

		real(8) g, height, dt, partial(1:2)
		integer(4) face_idx, x, y, dim, i, j

		g = var_pr.g;  height = var_pr.height;  dim = var_pr.dim

		dt = 10000d0


		do face_idx = 1, 6
			do y = -dim, dim
				do x = -dim, dim

					partial(1) = partial_c_x(var_pr.h_height(x-1:x+1, y, face_idx), var_pr.distance_grid(:, y, x, face_idx))
					var.u_vel(x, y, face_idx) = var_pr.u_vel(x, y, face_idx) - dt*g*partial(1)

					partial(1) = partial_c_y(var_pr.h_height(x, y-1:y+1, face_idx), var_pr.distance_grid(:, y, x, face_idx))
					var.v_vel(x, y, face_idx) = var_pr.v_vel(x, y, face_idx) - dt*g*partial(1)

					partial(1) = partial_c_x(var_pr.u_vel(x-1:x+1, y, face_idx), var_pr.distance_grid(:, y, x, face_idx))
					partial(2) = partial_c_y(var_pr.v_vel(x, y-1:y+1, face_idx), var_pr.distance_grid(:, y, x, face_idx))
					var.h_height(x, y, face_idx) = var_pr.h_height(x, y, face_idx) - height*(partial(1) + partial(2))

! 					if(face_idx == 2 .or. face_idx == 4) then
! 						if(abs(var.h_height(x, y, face_idx)) > height+10) var.h_height(x, y, face_idx) = 0
! 					else
! 						if(abs(var.h_height(x, y, face_idx)) > 50) var.h_height(x, y, face_idx) = 0
! 					end if


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

		call f.equal(var_pr, var)
		call f.borders(var_pr, var_pr)
		call f.corner_zero(var_pr)


	end subroutine



	real(8) function partial_c_x(fun, param)
		real(8), intent(in) :: fun(1:3), param(1:4)
		partial_c_x = (fun(3) - fun(1))/(param(2) + param(4))
	end function

	real(8) function partial_c_y(fun, param)
		real(8), intent(in) :: fun(1:3), param(1:4)
		partial_c_y = (fun(3) - fun(1))/(param(1) + param(3))
	end function



end module
