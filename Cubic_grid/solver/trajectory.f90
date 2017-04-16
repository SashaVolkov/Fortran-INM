module trajectory

	use grid_var, Only: g_var

implicit none

	Private
	Public :: traj

	Type traj

		Real(8), Allocatable :: alpha_x(:, :), alpha_y(:, :)
		Integer(4) ns_xy(2), nf_xy(2)


		CONTAINS
		Procedure, Public :: SemiLagr => SemiLagr

	End Type


CONTAINS


	Subroutine SemiLagr(this, g, interp, f, f_prev, f_next, mass)

		Class(Lagr) :: this
		Class(func) :: f, f_prev, f_next
		Class(intp) :: interp

		Real(8) ::  delta, x_wave, u_cov, v_cov
		Integer(4) :: x, i, p


		do face = 1,6
			do y = ns_y, nf_y
				do x = ns_x, nf_x
					u_cov(x, y, face, 0) = (3d0/2d0)*f.u_cov(x, y, face) - (1d0/2d0)*f_prev.u_cov(x, y, face)
					v_cov(x, y, face, 0) = (3d0/2d0)*f.v_cov(x, y, face) - (1d0/2d0)*f_prev.v_cov(x, y, face)
				end do
			end do
		end do


		do face = 1,6
			do y = ns_y, nf_y
				do x = ns_x, nf_x
					alpha_x(x, y, face, 0) = u_cov(x, y, face)*dt
					alpha_y(x, y, face, 0) = v_cov(x, y, face)*dt

					do i = 1, this.iter_stpes+1

						x_wave = x - alpha_x(x, y, face, i-1)
						y_wave = y - alpha_y(x, y, face, i-1)

							call interp.Lintp(x_wave, this.velocity, this.mass_alpha(i), g)

					end do
				end do
			end do
		end do

	End Subroutine


end module