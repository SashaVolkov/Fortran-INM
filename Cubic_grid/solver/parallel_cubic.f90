module parallel_cubic

implicit none

	include"mpif.h"

	Private
	Public :: parallel

	Type parallel
		integer(4) Ydim_block, Xdim_block, Xsize, Ysize, block_x, block_y, dim
		integer(4) ns_xy(1:2), nf_xy(1:2), step, x_grey, y_grey
		integer(4) snd_xy_up(1:2), snd_xy_down(1:2), snd_xy_right(1:2), snd_xy_left(1:2)
		integer(4) Neighbour(1:6, 1:4), border(6, 4), Neighbours_face(6, 4)
		CONTAINS
			Procedure, Public :: init => parallel_init
			Procedure, Private :: grey_zone => grey_zone
			Procedure, Private :: Neighbourhood => Neighbourhood
	End Type


CONTAINS


	Subroutine parallel_init(this, dim, step, np, id)

		Class(parallel) :: this

		integer, Intent(In) :: dim, step, np, id
		integer k, i, j, rc, p, ier, face
		integer dims(2)

		k=1; i = 0; j = 0; p = 1
		this.step = step
		this.dim = dim


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


		this.Xsize = (2*dim+1)/this.Xdim_block
		this.Ysize = (2*dim+1)/this.Ydim_block


		this.ns_xy(1) = 1 + this.Xsize*(this.block_x - 1); this.nf_xy(1) = this.Xsize*this.block_x
		this.ns_xy(2) = 1 + this.Ysize*(this.block_y - 1); this.nf_xy(2) = this.Ysize*this.block_y
! 		first_xy() = ns_xy(1) - step; last_xy = nf_xy(1) + step
! 		first_xy = ns_xy(2) - step; last_xy = nf_xy(2) + step


		call this.Neighbourhood(id)

		this.Xsize = 1 + this.nf_xy(1) - this.ns_xy(1)
		this.Ysize = 1 + this.nf_xy(2) - this.ns_xy(2)

		call this.grey_zone()

	End Subroutine



	Subroutine grey_zone(this)

		Class(parallel) :: this
		integer(4) x, y, mp_dp, ier


		x=this.Xsize + 2*this.step
		y=this.Ysize + 2*this.step

		mp_dp = MPI_DOUBLE

		call MPI_TYPE_VECTOR(this.Xsize, this.step, y, mp_dp, this.x_grey, ier)
		call MPI_TYPE_VECTOR(this.step, this.Ysize, x, mp_dp, this.y_grey, ier)
		call MPI_TYPE_COMMIT(this.x_grey, ier)
		call MPI_TYPE_COMMIT(this.y_grey, ier)

		this.snd_xy_left(1) = this.ns_xy(1) - this.step; this.snd_xy_left(2) = this.ns_xy(2);
		this.snd_xy_right(1) = this.nf_xy(1); this.snd_xy_right(2) = this.ns_xy(2);

		this.snd_xy_down(1) = this.ns_xy(1); this.snd_xy_down(2) = this.ns_xy(2) - this.step;
		this.snd_xy_up(1) = this.ns_xy(1); this.snd_xy_right(2) = this.nf_xy(2);

	End Subroutine



	Subroutine Neighbourhood(this, id)

	Class(parallel) :: this
	integer(4), Intent(In) :: id
	integer(4) face

	! 			Neighbourhood
	this.border(:, :) = 0 ! if 0 element in center, if "1": +pi/2 rotation, if "2": +pi, if "-1": -pi/2

	do face = 1, 6
		if (this.block_x < this.Xdim_block) then  ! right
			this.Neighbour(face, 2) = this.Ydim_block + id
			this.Neighbours_face(face, 2) = face
		else
			this.nf_xy(1) = 2*this.dim
			this.Neighbour(face, 2) = this.block_y - 1
			this.Neighbours_face(face, 2) = face+1
			if ( face == 1 ) then
				this.Neighbour(1, 2) = this.Ydim_block*(this.block_y - 1) + this.block_x - 1
				this.border(face, 2) = 1
				this.Neighbours_face(face, 2) = 3
			else if ( face == 6 ) then
				this.Neighbour(6, 2) = (this.Ydim_block - this.block_y)*this.Xdim_block
				this.border(face, 2) = -1
				this.Neighbours_face(face, 2) = 3
			end if
		end if


		if (this.block_x > 1) then  ! left
			this.Neighbour(face, 4) = id - this.Ydim_block
			this.Neighbours_face(face, 4) = face
		else
			this.ns_xy(1) = 1
			this.Neighbour(face, 4) = this.Ydim_block*(this.Xdim_block - 1) + this.block_y - 1
						if ( face == 1 ) then
							this.Neighbour(1, 4) = (this.Ydim_block - this.block_y+1)*this.Xdim_block - 1
							this.border(face, 4) = -1
						else if ( face == 6 ) then
							this.Neighbour(6, 4) = (this.block_y - 1)*this.Ydim_block
							this.border(face, 4) = 1
						end if
		end if


		if (this.block_y < this.Ydim_block) then  ! bottom
			this.Neighbour(face, 3) = 1 + id
			this.Neighbours_face(face, 3) = face
		else
			this.nf_xy(2) = 2*this.dim
			this.Neighbour(face, 3) = this.Ydim_block*(this.block_x-1)
			if ( face == 3 ) then
				this.Neighbour(3, 3) = this.Ydim_block*(this.Xdim_block - 1) + this.block_x - 1
				this.border(face, 3) = -1
			else if ( face == 5 ) then
				this.Neighbour(5, 3) = this.Xdim_block - this.block_x
				this.border(face, 3) = 1
			else if ( face == 4 ) then
				this.Neighbour(4, 3) = this.Ydim_block*(this.Xdim_block - this.block_x) + this.block_y -1
				this.border(face, 3) = 2
			else if ( face == 1 ) then
				this.border(face, 3) = 2
			end if
		end if


		if (this.block_y > 1) then  ! top
			this.Neighbour(face, 1) = id - 1
			this.Neighbours_face(face, 1) = face
		else
			this.ns_xy(2) = 1
			this.Neighbour(face, 1) = id + this.Ydim_block - 1
			if ( face == 3 ) then
				this.Neighbour(3, 1) = this.Ydim_block*this.Xdim_block - this.block_x
				this.border(face, 1) = 1
			else if ( face == 5 ) then
				this.Neighbour(5, 1) = this.block_x - 1; this.border(face, 1) = -1
			else if ( face == 4 ) then
				this.Neighbour(4, 1) = this.Ydim_block*(this.Xdim_block - this.block_x)
				this.border(face, 1) = 2
			else if ( face == 6 ) then
				this.border(face, 1) = 2
			end if
		end if
	end do


	End Subroutine



end module