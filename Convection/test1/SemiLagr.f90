Module SemiLagr

	Use modnet
	Use modfunc
	Use Interp
	Use MPI

IMPLICIT NONE

	Public :: Lagr

	Type Lagr

		Integer(4) lcase, eqvtype, iter_stpes


		Real(8), Allocatable :: velocity(:)
		Real(8), Allocatable :: xpoints(:)
		Real(8), Allocatable :: mass_alpha(:)

		CONTAINS

			Procedure ::  SemiLagr=> Lagr_SemiLagr
			Procedure ::  Interpol=> Lagr_Interpol

			Procedure ::  init=> Lagr_init
			Procedure ::  deinit=> Lagr_deinit


	End Type

	CONTAINS

	Subroutine Lagr_init(this, lcase, iter_stpes, eqvtype, g)

		Class(grid) :: g
		Class(Lagr) :: this
		Integer(4), Intent(In) :: iter_stpes, lcase, eqvtype
		Integer(4) x, y

		this.lcase = lcase
		this.iter_stpes = iter_stpes
		this.eqvtype = eqvtype

		Allocate(this.velocity_x(1:g.StepsY, g.ns-g.bstep:g.nf+g.fstep))
		Allocate(this.velocity_y(1:g.StepsY, g.ns-g.bstep:g.nf+g.fstep))
		Allocate(this.xpoints(1:g.StepsY, g.ns-g.bstep:g.nf+g.fstep))
		Allocate(this.mass_alpha(0:iter_stpes))


		do y = 1,g.StepsY
			do x = g.ns, g.nf
				this.velocity_x(y,x) = 0
				this.velocity_y(y,x) = 0
			end do
		end do


		do x = 0, iter_stpes
			this.mass_alpha(x) = 0
		end do


	End Subroutine


	Subroutine Lagr_SemiLagr(this, prev_mass, mass, g, interp)

		Class(Lagr) :: this
		Class(grid) :: g
		Class(intp) :: interp

		Real(8) ::  delta, x_wave
		Integer(4) :: x, i, p
		Real(8), Intent(inout) :: prev_mass(1:g.StepsY, g.ns - g.bstep : g.nf + g.fstep)
		Real(8), Intent(inout) :: mass(1:g.StepsY, g.ns - g.bstep : g.nf + g.fstep)


		do y = 1,g.StepsY
			do x = g.ns - g.bstep, g.nf + g.fstep
				if ( this.eqvtype == 2 ) then
					this.velocity_x(y,x) = prev_mass(y,x)
					this.velocity_y(y,x) = prev_mass(y,x)
				elseif ( this.eqvtype == 1 ) then
					this.velocity_x(y,x) = g.velocity_x
					this.velocity_y(y,x) = g.velocity_y
				end if
			end do
		end do





		do x = g.ns, g.nf
			this.mass_alpha(0) = this.velocity(x)

			do i = 1, this.iter_stpes+1

				x_wave = x - this.mass_alpha(i-1)
				delta = x_wave - int(x_wave)

				if ( i == this.iter_stpes+1 ) then
! 					call this.Interpol(delta, x_wave, mass(x), prev_mass, g)
					call interp.Lintp(x_wave, mass(x), prev_mass, g)
				else
! 					call this.Interpol(delta, x_wave, this.mass_alpha(i), this.velocity, g)
					call interp.Lintp(x_wave, this.mass_alpha(i), this.velocity, g)
				end if

			end do
		end do

	End Subroutine

	Subroutine Lagr_Interpol(this, delta, x_wave, output, mass_in, g)

		Class(Lagr) :: this
		Class(grid) :: g
		Real(8), Intent(In) :: x_wave
		Real(8), Intent(out) :: output
		Real(8), Intent(in) :: delta
		Real(8), Intent(In) :: mass_in(g.ns - g.bstep : g.nf + g.fstep)
		Integer(4) :: ix_wave
		ix_wave = int(x_wave)

		if ( this.lcase == 2 ) then
			output = (1.0 - delta)*mass_in(ix_wave) + delta*mass_in(ix_wave-1)
		elseif ( this.lcase == 3 ) then
			output = (delta/2)*(1+delta)*mass_in(ix_wave-1) + (1- delta**2)*mass_in(ix_wave) - (delta/2)*(1-delta)*mass_in(ix_wave+1)
		elseif ( this.lcase == 4 ) then
			output = -(delta/6)*(1-delta**2)*mass_in(ix_wave-2) + (delta/2)*(1+delta)*(2-delta)*mass_in(ix_wave-1) + ((1- delta**2)*(2-delta)*mass_in(ix_wave))/2 - (delta/6)*(1-delta)*(2-delta)*mass_in(ix_wave+1)
		end if

	End Subroutine

		Subroutine Lagr_deinit(this)
			Class(Lagr) :: this

			if (Allocated(this.velocity)) Deallocate(this.velocity)
			if (Allocated(this.mass_alpha)) Deallocate(this.mass_alpha)

		End Subroutine




End Module