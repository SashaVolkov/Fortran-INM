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

		Subroutine grid_init(this, LX, LY, LT, StepsX, StepsY, StepsT, bstep, fstep, velocity_x, velocity_y, np, id)

			Class(grid) :: this

			integer, Intent(In) :: StepsX, StepsY, StepsT, np, id, bstep, fstep
			real(8), Intent(In) :: LX, LY, LT, velocity_x, velocity_y


			this.LX = LX
			this.LY = LY
			this.LT = LT
			this.StepsX = StepsX
			this.StepsY = StepsY
			this.StepsT = StepsT

			this.bstep = bstep
			this.fstep = fstep
			this.velocity_x = velocity_x
			this.velocity_y = velocity_y

			this.np = np
			this.id = id


			this.dx = this.LX/this.StepsX
			this.dy = this.LY/this.StepsY
			this.dt = this.LT/this.StepsT
			this.Xsize = this.StepsX/np

			this.ns_x = id*this.Xsize+1
			this.nf_x = (id+1)*this.Xsize
			if ( id+1 == np ) then
				this.nf_x = this.StepsX
			end if
			
			this.ns_y = 1; this.nf_y = StepsY;


			this.Gamma_x = (this.velocity_x*this.dt)/this.dx
			this.Gamma_y = (this.velocity_y*this.dt)/this.dy

! 			print *, this.dy, this.dx

		End Subroutine


	End Module