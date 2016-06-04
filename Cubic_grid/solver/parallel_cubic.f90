module parallel_cubic

	! use grid_var, Only: g_var

implicit none

	include"mpif.h"

	Private
	Public :: parallel

	Type parallel
		integer(4) Neighbour(1:6, 1:4), Ydim_block, Xdim_block, Xsize, Ysize, block_x, block_y
		integer(4) ns_xy(1:2), nf_xy(1:2), step
		CONTAINS
			Procedure :: init => parallel_init
	End Type


CONTAINS


	Subroutine parallel_init(this, dim, step, np, id)

		Class(parallel) :: this

		integer, Intent(In) :: dim, step, np, id
! 		integer, Intent(Out) :: ns_xy(1:2), nf_xy(1:2)
		integer k, i, j, rc, p, ier, face
		integer dims(2)

		k=1; i = 0; j = 0; p = 1
		this.step = step


		do while ( np > k )
			k = k*2
			if (i > j) then
				j=j+1
			else
				i=i+1
			end if
		end do

		call MPI_DIMS_CREATE(np, 2, dims, ier)

		this.Ydim_block = dims(1); this.Xdim_block = dims(2)

			if ( id == 0 ) print *, this.Xdim_block, this.Ydim_block
			this.block_y = mod((id), this.Ydim_block)+1

			do while ( id+1 > p*this.Ydim_block )
				p = p + 1
			end do
			this.block_x = p
! 			print *,"table", this.block_x, this.block_y, id+1
! 		end if


		this.Xsize = (2*dim+1)/this.Xdim_block
		this.Ysize = (2*dim+1)/this.Ydim_block



		this.ns_xy(1) = 1 + this.Xsize*(this.block_x - 1) - dim; this.nf_xy(1) = this.Xsize*this.block_x - dim
		this.ns_xy(2) = 1 + this.Ysize*(this.block_y - 1) - dim; this.nf_xy(2) = this.Ysize*this.block_y - dim
! 		first_xy() = ns_xy(1) - step; last_xy = nf_xy(1) + step
! 		first_xy = ns_xy(2) - step; last_xy = nf_xy(2) + step


! 			Neighbourhood
	do face = 1, 6
		if (this.block_x < this.Xdim_block) then  ! right
			this.Neighbour(face, 2) = this.Ydim_block + id
		else
			this.nf_xy(1) = dim
			this.Neighbour(face, 2) = this.block_y - 1
			if ( face == 1 ) then
				this.Neighbour(1, 2) = this.Ydim_block*(this.block_y - 1) + this.block_x - 1
			else if ( face == 6 ) then
				this.Neighbour(6, 2) = (this.Ydim_block - this.block_y)*this.Xdim_block
			end if
		end if

		if (this.block_x > 1) then  ! left
			this.Neighbour(face, 4) = id - this.Ydim_block
		else
			this.ns_xy(1) = -dim
			this.Neighbour(face, 4) = this.Ydim_block*(this.Xdim_block - 1) + this.block_y - 1
						if ( face == 1 ) then
							this.Neighbour(1, 4) = (this.Ydim_block - this.block_y+1)*this.Xdim_block - 1
						else if ( face == 6 ) then
							this.Neighbour(6, 4) = (this.block_y - 1)*this.Ydim_block
						end if
		end if


		if (this.block_y < this.Ydim_block) then  ! bottom
			this.Neighbour(face, 3) = 1 + id
		else
			this.nf_xy(2) = dim
			this.Neighbour(face, 3) = this.Ydim_block*(this.block_x-1)
			if ( face == 3 ) then
				this.Neighbour(3, 3) = this.Ydim_block*(this.Xdim_block - 1) + this.block_x - 1
			else if ( face == 5 ) then
				this.Neighbour(5, 3) = this.Xdim_block - this.block_x
			else if ( face == 4 ) then
				this.Neighbour(4, 3) =  this.Ydim_block*(this.Xdim_block - this.block_x) + this.block_y -1
			end if
		end if

		if (this.block_y > 1) then  ! top
			this.Neighbour(face, 1) = id - 1
		else
			this.ns_xy(2) = -dim
			this.Neighbour(face, 1) = id + this.Ydim_block - 1
			if ( face == 3 ) then
				this.Neighbour(3, 1) = this.Ydim_block*this.Xdim_block - this.block_x
			else if ( face == 5 ) then
				this.Neighbour(5, 1) = this.block_x - 1
			else if ( face == 4 ) then
				this.Neighbour(4, 1) = this.Ydim_block*(this.Xdim_block - this.block_x)
			end if
		end if
	end do


	this.Xsize = 1 + this.nf_xy(1) - this.ns_xy(1)
	this.Ysize = 1 + this.nf_xy(2) - this.ns_xy(2)

	! if ( id == 3) then
	! 	print *, this.Neighbour(4, :)
	! 	print *, this.ns_xy, this.nf_xy
	! 	print *, this.Xsize, this.Ysize
	! end if


	End Subroutine



end module