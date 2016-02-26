Module modnet

	IMPLICIT NONE

	Public :: grid

	Type grid

		real(8) LX, LY, LT, velocity_x, velocity_y, dx, dy, dt, Gamma_x, Gamma_y

		integer StepsX, StepsY, StepsT, Xsize, ns_x,nf_x, Ysize, ns_y,nf_y, np, id, bstep, fstep, block_x, block_y, Xdim_block, Ydim_block
		integer Neighb_up, Neighb_down, Neighb_left, Neighb_right


		CONTAINS
			Procedure :: init => grid_init

	End Type

	CONTAINS

	Subroutine grid_init(this, LX, LY, LT, StepsX, StepsY, StepsT, bstep, fstep, np, id)

		Class(grid) :: this

		integer, Intent(In) :: StepsX, StepsY, StepsT, np, id, bstep, fstep
		real(8), Intent(In) :: LX, LY, LT
		integer k, i, j, rc, p

		k=1; i = 0; j = 0; p = 1

		this.LX = LX; this.LY = LY; this.LT = LT

		this.StepsX = StepsX; this.StepsY = StepsY; this.StepsT = StepsT

		this.bstep = bstep; this.fstep = fstep

		this.np = np; this.id = id

		this.dx = this.LX/this.StepsX
		this.dy = this.LY/this.StepsY
		this.dt = this.LT/this.StepsT

		do while ( np > k )
			k = k*2
			if (i > j) then
				j=j+1
			else
				i=i+1
			end if
		end do

		this.Ydim_block = 2**(j); this.Xdim_block = 2**(i)

		if ( k > np ) then
			if ( id == 0 ) print *, "Use number of process 2^n", this.Xdim_block, this.Ydim_block
			call MPI_FINALIZE(rc)
			STOP
		else
			if ( id == 0 ) print *, this.Xdim_block, this.Ydim_block
			this.block_x = mod((id), this.Ydim_block)+1

			do while ( id+1 > p*this.Ydim_block )
				p = p + 1
			end do
			this.block_y = p
! 			print *,"table", this.block_x, this.block_y, id+1
		end if

		if (np == 1) then
			this.Xsize = this.StepsX
			this.Ysize = this.StepsY
		else
			this.Xsize = this.StepsX/this.Xdim_block
			this.Ysize = this.StepsY/this.Ydim_block
		end if
		
		this.ns_x = 1; this.nf_x = this.StepsX
		this.ns_y = 1; this.nf_y = StepsY


! 			Neighbourhood

	if (this.block_y < this.Xdim_block) then
		this.Neighb_right = this.Ydim_block + id
	else
		this.Neighb_right = -1
	end if

	if (this.block_y > 1) then
		this.Neighb_left = id - this.Ydim_block
	else
		this.Neighb_left = -1
	end if


	if (this.block_x < this.Ydim_block) then
		this.Neighb_down = 1 + id
	else
		this.Neighb_down = -1
	end if

	if (this.block_x > 1) then
		this.Neighb_up = id - 1
	else
		this.Neighb_up = -1
	end if


	if ( id == 5 ) print *,"id:",id,"right:",this.Neighb_right,"left:",this.Neighb_left,"down:",this.Neighb_down,"up:",this.Neighb_up


	End Subroutine


End Module