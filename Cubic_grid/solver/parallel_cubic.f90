module parallel_cubic

implicit none

	include"mpif.h"

	Private
	Public :: parallel

	Type parallel
		integer(4) Ydim_block, Xdim_block, Xsize, Ysize, block_x, block_y, dim, first_x, first_y, last_x, last_y
		integer(4) ns_xy(1:2), nf_xy(1:2), step, up, right, left, down, halo(6, 4), y_vec_type, rot(6, 4)
		integer(4) snd_xy(6, 4, 2), snd_xy_180(6, 4, 2), snd_xy_90(6, 4, 2), snd_xy_m90(6, 4, 2)
		integer(4) rcv_xy(6, 4, 2), rcv_xy_90(6, 4, 2), rcv_xy_m90(6, 4, 2)
		integer(4) Neighbour_id(1:6, 1:4), border(6, 4), Neighbours_face(6, 4), id
		CONTAINS
			Procedure, Public :: init => parallel_init
			Procedure, Private :: halo_zone => halo_zone
			Procedure, Private :: Neighbourhood => Neighbourhood
	End Type


CONTAINS



! For example, in C, a two-dimensional array with three rows and four columns will be stored in memory in the following sequence:
! (1,1),(1,2),(1,3),(1,4),(2,1),(2,2),(2,3),(2,4),(3,1),(3,2),(3,3),(3,4)

! However, FORTRAN will store the same array in the following sequence:
! (1,1),(2,1),(3,1),(1,2),(2,2),(3,2),(1,3),(2,3),(3,3),(1,4),(2,4),(3,4)



	Subroutine parallel_init(this, dim, step, np, id)

		Class(parallel) :: this

		integer, Intent(In) :: dim, step, np, id
		integer k, i, j, rc, p, ier, face
		integer dims(2)

		k=1; i = 0; j = 0; p = 0
		this.step = step
		this.dim = dim
		this.id = id

		this.up = 1; this.right=2; this.down=3; this.left=4


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
		this.block_y = mod((id), this.Ydim_block)

		do while ( id > (p + 1)*this.Ydim_block - 1)
			p = p + 1
		end do
		this.block_x = p


		this.Xsize = 2*dim/this.Xdim_block
		this.Ysize = 2*dim/this.Ydim_block


		this.ns_xy(1) = 1 + this.Xsize*this.block_x  ; this.nf_xy(1) = this.Xsize*(this.block_x + 1)
		this.ns_xy(2) = 1 + this.Ysize*this.block_y  ; this.nf_xy(2) = this.Ysize*(this.block_y + 1)

		this.first_x = this.ns_xy(1) - step;  this.first_y = this.ns_xy(2) - step
		this.last_x = this.nf_xy(1) + step;  this.last_y = this.nf_xy(2) + step


		call this.Neighbourhood(id)

		this.Xsize = 1 + this.nf_xy(1) - this.ns_xy(1)
		this.Ysize = 1 + this.nf_xy(2) - this.ns_xy(2)

		do face = 1, 6
			call this.halo_zone(face)
		end do

	End Subroutine



	Subroutine halo_zone(this, face)

		Class(parallel) :: this
		integer(4), intent(in) :: face
		integer(4) x, y, mp_dp, ier, i, k, n
		integer(4) up, right, left, down, displ(1:this.step*this.Ysize), blocklen(1:this.step*this.Ysize)
		
		up = 1; right=2; down=3; left=4
		x=this.Xsize + 2*this.step
		mp_dp = MPI_DOUBLE_PRECISION
		n = this.step*this.Ysize

		displ(1) = 0
		blocklen(1) = 1
		do k = 2, this.Ysize
			displ(k) = displ(k - 1) + x
			blocklen(k) = 1
		end do

		if(this.step > 1) then
			do i = this.Ysize+1, n
				displ(i) = displ(i - this.Ysize) + 1
				blocklen(i) = 1
			end do
		end if


		this.rcv_xy(face, left, 1) = this.ns_xy(1) - this.step
		this.rcv_xy(face, left, 2) = this.ns_xy(2);
		this.rcv_xy(face, right, 1) = this.nf_xy(1) + 1
		this.rcv_xy(face, right, 2) = this.ns_xy(2);

		this.rcv_xy(face, down, 1) = this.ns_xy(1)
		this.rcv_xy(face, down, 2) = this.ns_xy(2) - this.step;
		this.rcv_xy(face, up, 1) = this.ns_xy(1)
		this.rcv_xy(face, up, 2) = this.nf_xy(2) + 1;


		this.snd_xy(face, left, :) = this.ns_xy(:)
		this.snd_xy(face, right, 1) = this.nf_xy(1) - this.step + 1
		this.snd_xy(face, right, 2) = this.ns_xy(2)

		this.snd_xy(face, down, :) = this.ns_xy(:)
		this.snd_xy(face, up, 1) = this.ns_xy(1)
		this.snd_xy(face, up, 2) = this.nf_xy(2) - this.step + 1;


		call MPI_TYPE_VECTOR(this.step, this.Xsize, x, mp_dp, this.halo(face, up), ier)
		call MPI_TYPE_VECTOR(this.step, this.Xsize, x, mp_dp, this.halo(face, down), ier)
		call MPI_TYPE_INDEXED(n, blocklen, displ, mp_dp, this.halo(face, right), ier)
		call MPI_TYPE_INDEXED(n, blocklen, displ, mp_dp, this.halo(face, left), ier)

		do i = 1, 4
			call MPI_TYPE_COMMIT(this.halo(face, i), ier)
		end do


	End Subroutine



	Subroutine Neighbourhood(this, id)

	Class(parallel) :: this
	integer(4), Intent(In) :: id
	integer(4) face, up, right, left, down, i
	
	up = 1; right=2; down=3; left=4

	! 			Neighbourhood
	this.border(:, :) = 0 ! if 0 element in center, if "1": +pi/2 rotation, if "2": +pi, if "-1": -pi/2
	this.Neighbour_id(:, 1) = id + 1
	this.Neighbour_id(:, 2) = this.Ydim_block + id
	this.Neighbour_id(:, 3) = id - 1
	this.Neighbour_id(:, 4) = id - this.Ydim_block
	this.rot(:, :) = 0

	do face = 1, 6
			this.Neighbours_face(face, :) = face

		if (this.block_y == this.Ydim_block - 1) then  ! up

			this.Neighbour_id(face, up) = this.block_x*this.Xdim_block
			this.Neighbours_face(face, up) = 6

			if ( face == 3 ) then
				this.Neighbour_id(face, up) = this.block_y*this.Xdim_block + this.block_x
				this.border(face, up) = -1
				this.rot(face, up) = 1

			else if ( face == 5 ) then
				this.Neighbour_id(face, up) = this.Xdim_block - this.block_x - 1
				this.border(face, up) = 1
				this.rot(face, up) = -1

			else if ( face == 4 ) then
				this.Neighbour_id(face, up) = this.Ydim_block*this.Xdim_block - 1 - this.block_x*this.Xdim_block
				this.border(face, up) = 2

			else if ( face == 6 ) then
				this.Neighbour_id(face, up) = this.Ydim_block*this.Xdim_block - 1 - this.block_x*this.Xdim_block
				this.Neighbours_face(face, up) = 4
				this.border(face, up) = 2

			else if ( face == 1 ) then
				this.Neighbours_face(face, up) = 2

			end if
		end if



		if (this.block_x == this.Xdim_block - 1) then  ! right

			this.Neighbour_id(face, right) = this.block_y
			this.Neighbours_face(face, right) = face + 1

			if ( face == 1 ) then
				this.Neighbour_id(face, right) = (this.Ydim_block - this.block_y - 1)*this.Xdim_block
				this.border(face, right) = 1
				this.rot(face, right) = 1
				this.Neighbours_face(face, right) = 3

			else if ( face == 6 ) then
				this.Neighbour_id(face, right) = this.block_y*this.Xdim_block + this.block_x
				this.border(face, right) = -1
				this.rot(face, right) = -1
				this.Neighbours_face(face, right) = 3

			else if(face == 5) then
				this.Neighbours_face(face, right) = 2
			end if
		end if


		if (this.block_y == 0) then  ! down

			this.Neighbour_id(face, down) = this.block_x*this.Xdim_block + this.Ydim_block - 1
			this.Neighbours_face(face, down) = 1

			if ( face == 3 ) then
				this.Neighbour_id(face, down) = this.Ydim_block*this.Xdim_block - this.block_x - 1
				this.border(face, down) = 1
				this.rot(face, down) = -1

			else if ( face == 5 ) then
				this.Neighbour_id(face, down) = this.block_x
				this.border(face, down) = -1
				this.rot(face, down) = 1

			else if ( face == 4 ) then
				this.Neighbour_id(face, down) = this.Ydim_block*(this.Xdim_block - 1) - this.block_x*this.Xdim_block
				this.border(face, down) = 2

			else if ( face == 1 ) then
				this.Neighbour_id(face, down) = this.Ydim_block*(this.Xdim_block - 1) - this.block_x*this.Xdim_block
				this.border(face, down) = 2
				this.Neighbours_face(face, down) = 4

			else if ( face == 6 ) then
				this.Neighbours_face(face, down) = 2
			end if
		end if


		if (this.block_x == 0) then  ! left

			this.Neighbour_id(face, left) = this.block_y + this.Xdim_block*(this.Xdim_block - 1)
			this.Neighbours_face(face, left) = face - 1

			if ( face == 1 ) then
				this.Neighbour_id(face, left) = this.Xdim_block*this.block_y
				this.border(face, left) = -1
				this.rot(face, left) = -1
				this.Neighbours_face(face, left) = 5

			else if ( face == 6 ) then
				this.Neighbour_id(face, left) = this.Xdim_block*this.Ydim_block - 1 - this.block_y*this.Xdim_block
				this.border(face, left) = 1
				this.rot(face, left) = 1
				this.Neighbours_face(face, left) = 5
			else if ( face == 2 ) then
				this.Neighbours_face(face, left) = 5
			end if
		end if
	end do


	End Subroutine



end module