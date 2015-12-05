	Module modnet

		IMPLICIT NONE

		Public :: grid

		Type grid

			real(8) LX, LT, velocity, dx, dt, Gamma

			integer StepsX, StepsT, Xsize, ns,nf, np, id, bstep, fstep

			CONTAINS
				Procedure :: init => grid_init

		End Type

		CONTAINS

		Subroutine grid_init(this, LX, LT, StepsX, StepsT, bstep, fstep, velocity, np, id)

			Class(grid) :: this

			integer, Intent(In) :: StepsX, StepsT, np, id, bstep, fstep
			real(8), Intent(In) :: LX, LT, velocity


			this.LX = LX
			this.LT = LT
			this.StepsX = StepsX
			this.StepsT = StepsT

			this.bstep = bstep
			this.fstep = fstep
			this.velocity = velocity

			this.np = np
			this.id = id


			this.dx = this.LX/this.StepsX
			this.dt = this.LT/this.StepsT
			this.Xsize = this.StepsX/np

			this.ns = id*this.Xsize+1
			this.nf = (id+1)*this.Xsize
			if ( id+1 == np ) then
				this.nf = this.StepsX
			end if
			
			this.Gamma = (this.Velocity*this.dt)/this.dx

		End Subroutine


	End Module