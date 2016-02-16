	Module modnet

		IMPLICIT NONE

		Public :: grid

		Type grid

			real(8) LX, LY, LT, velocity_x, velocity_y, dx, dy, dt, Gamma_x, Gamma_y

			integer StepsX, StepsY, StepsT, Xsize, ns_x,nf_x, Ysize, ns_y,nf_y, np, id, bstep, fstep

			CONTAINS
				Procedure :: init => grid_init

		End Type

		CONTAINS

		Subroutine grid_init(this, LX, LY, LT, StepsX, StepsY, StepsT, bstep, fstep, np, id)

			Class(grid) :: this

			integer, Intent(In) :: StepsX, StepsY, StepsT, np, id, bstep, fstep
			real(8), Intent(In) :: LX, LY, LT

			this.LX = LX
			this.LY = LY
			this.LT = LT
			this.StepsX = StepsX
			this.StepsY = StepsY
			this.StepsT = StepsT

			this.bstep = bstep
			this.fstep = fstep


			this.np = np
			this.id = id


			this.dx = this.LX/this.StepsX
			this.dy = this.LY/this.StepsY
			this.dt = this.LT/this.StepsT
			this.Xsize = this.StepsX/np

			this.ns_x = 1; this.nf_x = this.Xsize
			
			this.ns_y = 1; this.nf_y = StepsY;



		End Subroutine


	End Module