module parallel_cubic

	use mpi

implicit none


	Private
	Public :: parallel

	Type parallel
		Integer(4) Ydim_block, Xdim_block, Xsize, Ysize, block_x, block_y, dim, first_x, first_y, last_x, last_y
		Integer(4) ns_xy(1:2), nf_xy(1:2), step, up, right, left, down, halo(6, 4), halo_corn(6, 4), rot(6, 0:4)
		Integer(4) snd_xy(6, 4, 2), rcv_xy(6, 4, 2), corn_snd_xy(6, 4, 2), corn_rcv_xy(6, 4, 2)
		Integer(4) Neighbour_id(6, 4), border(6, 0:4), Neighbours_face(6, 4), id, Neighb_dir(6,4)
		Integer(4) Neighbour_corn_id(6, 4), Neighb_corn_dir(6,4), Neighbours_corn_face(6, 4), My_dir(6, 4)
		CONTAINS
			Procedure, Public :: init => parallel_init
			Procedure, Private :: halo_zone => halo_zone
			Procedure, Private :: Neighbourhood => Neighbourhood
			Procedure, Private :: Displacement => Displacement
			Procedure, Private :: Displacement_corn => Displacement_corn
	End Type


CONTAINS



! For example, in C, a two-dimensional array with three rows and four columns will be stored in memory in the following sequence:
! (1,1),(1,2),(1,3),(1,4),(2,1),(2,2),(2,3),(2,4),(3,1),(3,2),(3,3),(3,4)

! However, FORTRAN will store the same array in the following sequence:
! (1,1),(2,1),(3,1),(1,2),(2,2),(3,2),(1,3),(2,3),(3,3),(1,4),(2,4),(3,4)

! ______________________
! |   |             |   |
! | D |      up     | A |
! |___|_____________|___|
! |   |             | r |
! | l |             | i |
! | e |             | g |
! | f |             | h |
! | t |             | t |
! |___|_____________|___|
! |   |             |   |
! | C |     down    | B |
! |___|_____________|___|


Subroutine parallel_init(this, dim, step, np, id)

	Class(parallel) :: this

	Integer, Intent(In) :: dim, step, np, id
	Integer k, i, j, rc, p, ier, face
	Integer dims(2)

! 	call MPI_Comm_rank(MPI_COMM_WORLD,id,ier)
! 	call MPI_Comm_size(MPI_COMM_WORLD,np,ier)

	k=1; i = 0; j = 0; p = 0
	this.step = step
	! if(step == 1) this.step = 2
	this.dim = dim
	this.id = id

	this.up = 1; this.right=2; this.down=3; this.left=4


	! call MPI_DIMS_CREATE(np, 2, dims, ier)
	! print *, ier
	dims(1) = sqrt(dble(np))
	dims(2) = dims(1)

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

	this.first_x = this.ns_xy(1) - this.step;  this.first_y = this.ns_xy(2) - this.step
	this.last_x = this.nf_xy(1) + this.step;  this.last_y = this.nf_xy(2) + this.step


	call this.Neighbourhood(id)

	this.Xsize = 1 + this.nf_xy(1) - this.ns_xy(1)
	this.Ysize = 1 + this.nf_xy(2) - this.ns_xy(2)

	do face = 1, 6
		call this.halo_zone(face)
	end do

End Subroutine



Subroutine halo_zone(this, face)

	Class(parallel) :: this
	Integer(4), intent(in) :: face
	Integer(4) x, y, mp_dp, ier, i, j, k, n, A, B, C, D, displ_corn(1:this.step*this.step, 4)
	Integer(4) up, right, left, down, displ(1:this.step*this.Xsize), blocklen(1:this.step*this.Ysize)

	up = 1; right=2; down=3; left=4;  A=1;  B=2;  C=3;  D=4
	x=this.Xsize + 2*this.step
	mp_dp = MPI_DOUBLE_PRECISION
	n = this.step*this.Ysize


	this.rcv_xy(face, left, 1) = this.ns_xy(1) - this.step
	this.rcv_xy(face, left, 2) = this.ns_xy(2)
	this.rcv_xy(face, right, 1) = this.nf_xy(1) + 1
	this.rcv_xy(face, right, 2) = this.ns_xy(2)

	this.rcv_xy(face, down, 1) = this.ns_xy(1)
	this.rcv_xy(face, down, 2) = this.ns_xy(2) - this.step
	this.rcv_xy(face, up, 1) = this.ns_xy(1)
	this.rcv_xy(face, up, 2) = this.nf_xy(2) + 1

	this.snd_xy(face, left, :) = this.ns_xy(:)
	this.snd_xy(face, right, 1) = this.nf_xy(1) - this.step + 1
	this.snd_xy(face, right, 2) = this.ns_xy(2)

	this.snd_xy(face, down, :) = this.ns_xy(:)
	this.snd_xy(face, up, 1) = this.ns_xy(1)
	this.snd_xy(face, up, 2) = this.nf_xy(2) - this.step + 1



	this.corn_rcv_xy(face, A, 1) = this.nf_xy(1)+1
	this.corn_rcv_xy(face, A, 2) = this.nf_xy(2)+1
	this.corn_rcv_xy(face, B, 1) = this.nf_xy(1) + 1
	this.corn_rcv_xy(face, B, 2) = this.first_y

	this.corn_rcv_xy(face, C, 1) = this.ns_xy(1) - this.step
	this.corn_rcv_xy(face, C, 2) = this.ns_xy(2) - this.step
	this.corn_rcv_xy(face, D, 1) = this.ns_xy(1) - this.step
	this.corn_rcv_xy(face, D, 2) = this.nf_xy(2) + 1

	this.corn_snd_xy(face, A, 1) = this.nf_xy(1) - this.step + 1
	this.corn_snd_xy(face, A, 2) = this.nf_xy(2) - this.step + 1
	this.corn_snd_xy(face, B, 1) = this.nf_xy(1) - this.step + 1
	this.corn_snd_xy(face, B, 2) = this.ns_xy(2)

	this.corn_snd_xy(face, C, :) = this.ns_xy(:)
	this.corn_snd_xy(face, D, 1) = this.ns_xy(1)
	this.corn_snd_xy(face, D, 2) = this.nf_xy(2) - this.step + 1


	do i = 1, n
		blocklen(i) = 1
	end do

	do i = 1, 4

		call this.Displacement(face, i, displ)
		call MPI_TYPE_INDEXED(n, blocklen, displ, mp_dp, this.halo(face, i), ier)
		call MPI_TYPE_COMMIT(this.halo(face, i), ier)
		! if (face == 4 .and. i == 1 .and. this.id == 3) print *, displ(:)

	end do

	do i = 1, 4
		call this.Displacement_corn(face, i, displ_corn)
		call MPI_TYPE_INDEXED(this.step**2, blocklen(1:this.step**2), displ_corn, mp_dp, this.halo_corn(face, i), ier)
		call MPI_TYPE_COMMIT(this.halo_corn(face, i), ier)
	end do

End Subroutine



Subroutine Neighbourhood(this, id)

	Class(parallel) :: this
	Integer(4), Intent(In) :: id
	Integer(4) face, up, right, left, down, A, B, C, D, i
	
	up = 1; right=2; down=3; left=4;  A=1;  B=2;  C=3;  D=4

	! 			Neighbourhood
	this.border(:, :) = 0 ! if "0" no rotation, if "1": +pi/2 rotation, if "2": +pi, if "-1": -pi/2
	this.Neighbour_id(:, up) = id + 1
	this.Neighbour_id(:, right) = this.Ydim_block + id
	this.Neighbour_id(:, down) = id - 1
	this.Neighbour_id(:, left) = id - this.Ydim_block

	this.Neighbour_corn_id(:, A) = this.Ydim_block + id + 1
	this.Neighbour_corn_id(:, B) = this.Ydim_block + id - 1
	this.Neighbour_corn_id(:, C) = id - this.Ydim_block - 1
	this.Neighbour_corn_id(:, D) = id - this.Ydim_block + 1

	this.rot(:, :) = 0
	this.Neighb_dir(:, up) = down
	this.Neighb_dir(:, down) = up
	this.Neighb_dir(:, left) = right
	this.Neighb_dir(:, right) = left

	this.Neighb_corn_dir(:, A) = C
	this.Neighb_corn_dir(:, B) = D
	this.Neighb_corn_dir(:, C) = A
	this.Neighb_corn_dir(:, D) = B

	this.My_dir(:, A) = up
	this.My_dir(:, B) = down
	this.My_dir(:, C) = down
	this.My_dir(:, D) = up


	do face = 1, 6
			this.Neighbours_face(face, :) = face
			this.Neighbours_corn_face(face, :) = face

		if (this.block_y == this.Ydim_block - 1) then  ! up

			this.Neighbour_id(face, up) = this.block_x*this.Ydim_block
			this.Neighbour_corn_id(face, A) = (this.block_x + 1)*this.Ydim_block
			this.Neighbour_corn_id(face, D) = (this.block_x - 1)*this.Ydim_block
			this.Neighbours_face(face, up) = 6
			this.Neighbours_corn_face(face, A) = 6
			this.Neighbours_corn_face(face, D) = 6
			this.My_dir(face, A) = up
			this.My_dir(face, D) = up

			if ( face == 3 ) then
				this.Neighbour_id(face, up) = this.block_y*this.Xdim_block + this.block_x
				this.Neighbour_corn_id(face, A) = this.block_y*this.Xdim_block + this.block_x + 1
				this.Neighbour_corn_id(face, D) = this.block_y*this.Xdim_block + this.block_x - 1
				this.border(face, up) = -1
				this.rot(face, up) = 1
				this.Neighb_dir(face, up) = right
				this.Neighb_corn_dir(face, A) = B
				this.Neighb_corn_dir(face, D) = A

			else if ( face == 5 ) then
				this.Neighbour_id(face, up) = this.Xdim_block - this.block_x - 1
				this.Neighbour_corn_id(face, A) = this.Xdim_block - this.block_x - 2
				this.Neighbour_corn_id(face, D) = this.Xdim_block - this.block_x
				this.border(face, up) = 1
				this.rot(face, up) = -1
				this.Neighb_dir(face, up) = left
				this.Neighb_corn_dir(face, A) = D
				this.Neighb_corn_dir(face, D) = C

			else if ( face == 4 ) then
				this.Neighbour_id(face, up) = this.Ydim_block*this.Xdim_block - 1 - this.block_x*this.Xdim_block
				this.Neighbour_corn_id(face, A) = this.Ydim_block*this.Xdim_block - 1 - (this.block_x+1)*this.Xdim_block
				this.Neighbour_corn_id(face, D) = this.Ydim_block*this.Xdim_block - 1 - (this.block_x-1)*this.Xdim_block
				this.border(face, up) = 2
				this.Neighb_dir(face, up) = up
				this.Neighb_corn_dir(face, A) = A
				this.Neighb_corn_dir(face, D) = D

			else if ( face == 6 ) then
				this.Neighbour_id(face, up) = this.Ydim_block*this.Xdim_block - 1 - this.block_x*this.Xdim_block
				this.Neighbour_corn_id(face, A) = this.Ydim_block*this.Xdim_block - 1 - (this.block_x+1)*this.Xdim_block
				this.Neighbour_corn_id(face, D) = this.Ydim_block*this.Xdim_block - 1 - (this.block_x-1)*this.Xdim_block
				this.Neighbours_face(face, up) = 4
				this.Neighbours_corn_face(face, A) = 4
				this.Neighbours_corn_face(face, D) = 4
				this.rot(face, up) = 2
				this.border(face, up) = 2
				this.Neighb_dir(face, up) = up
				this.Neighb_corn_dir(face, A) = A
				this.Neighb_corn_dir(face, D) = D

			else if ( face == 1 ) then
				this.Neighbours_face(face, up) = 2
				this.Neighbours_corn_face(face, A) = 2
				this.Neighbours_corn_face(face, D) = 2

			end if
		end if



		if (this.block_x == this.Xdim_block - 1) then  ! right

			this.Neighbour_id(face, right) = this.block_y
			this.Neighbour_corn_id(face, A) = this.block_y + 1
			this.Neighbour_corn_id(face, B) = this.block_y - 1
			this.Neighbours_face(face, right) = face + 1
			this.Neighbours_corn_face(face, A) = face + 1
			this.Neighbours_corn_face(face, B) = face + 1
			this.My_dir(face, A) = right
			this.My_dir(face, B) = right

			if ( face == 1 ) then
				this.Neighbour_id(face, right) = (this.Ydim_block - this.block_y - 1)*this.Xdim_block
				this.Neighbour_corn_id(face, A) = (this.Ydim_block - this.block_y - 2)*this.Xdim_block
				this.Neighbour_corn_id(face, B) = (this.Ydim_block - this.block_y)*this.Xdim_block
				this.border(face, right) = 1
				this.rot(face, right) = 1
				this.Neighbours_face(face, right) = 3
				this.Neighbours_corn_face(face, A) = 3
				this.Neighbours_corn_face(face, B) = 3
				this.Neighb_dir(face, right) = down
				this.Neighb_corn_dir(face, A) = B
				this.Neighb_corn_dir(face, B) = C

			else if ( face == 6 ) then
				this.Neighbour_id(face, right) = this.block_y*this.Xdim_block + this.block_x
				this.Neighbour_corn_id(face, A) = (this.block_y + 1)*this.Xdim_block + this.block_x
				this.Neighbour_corn_id(face, B) = (this.block_y - 1)*this.Xdim_block + this.block_x
				this.border(face, right) = -1
				this.rot(face, right) = -1
				this.Neighbours_face(face, right) = 3
				this.Neighbours_corn_face(face, A) = 3
				this.Neighbours_corn_face(face, B) = 3
				this.Neighb_dir(face, right) = up
				this.Neighb_corn_dir(face, A) = D
				this.Neighb_corn_dir(face, B) = A

			else if(face == 5) then
				this.Neighbours_face(face, right) = 2
				this.Neighbours_corn_face(face, A) = 2
				this.Neighbours_corn_face(face, B) = 2
			end if
		end if


		if (this.block_y == 0) then  ! down

			this.Neighbour_id(face, down) = this.block_x*this.Xdim_block + this.Ydim_block - 1
			this.Neighbour_corn_id(face, B) = (this.block_x+1)*this.Xdim_block + this.Ydim_block - 1
			this.Neighbour_corn_id(face, C) = (this.block_x-1)*this.Xdim_block + this.Ydim_block - 1
			this.Neighbours_face(face, down) = 1
			this.Neighbours_corn_face(face, B) = 1
			this.Neighbours_corn_face(face, C) = 1
			this.My_dir(face, B) = down
			this.My_dir(face, C) = down

			if ( face == 3 ) then
				this.Neighbour_id(face, down) = this.Ydim_block*this.Xdim_block - this.block_x - 1
				this.Neighbour_corn_id(face, B) = this.Ydim_block*this.Xdim_block - this.block_x - 2
				this.Neighbour_corn_id(face, C) = this.Ydim_block*this.Xdim_block - this.block_x
				this.border(face, down) = 1
				this.rot(face, down) = -1
				this.Neighb_dir(face, down) = right
				this.Neighb_corn_dir(face, B) = A
				this.Neighb_corn_dir(face, C) = B

			else if ( face == 5 ) then
				this.Neighbour_id(face, down) = this.block_x
				this.Neighbour_corn_id(face, B) = this.block_x + 1
				this.Neighbour_corn_id(face, C) = this.block_x - 1
				this.border(face, down) = -1
				this.rot(face, down) = 1
				this.Neighb_dir(face, down) = left
				this.Neighb_corn_dir(face, B) = C
				this.Neighb_corn_dir(face, C) = D

			else if ( face == 4 ) then
				this.Neighbour_id(face, down) = this.Ydim_block*(this.Xdim_block - 1) - this.block_x*this.Xdim_block
				this.Neighbour_corn_id(face, B) = this.Ydim_block*(this.Xdim_block - 2) - this.block_x*this.Xdim_block
				this.Neighbour_corn_id(face, C) = this.Ydim_block*(this.Xdim_block) - this.block_x*this.Xdim_block
				this.border(face, down) = 2
				this.Neighb_dir(face, down) = down
				this.Neighb_corn_dir(face, B) = B
				this.Neighb_corn_dir(face, C) = C

			else if ( face == 1 ) then
				this.Neighbour_id(face, down) = this.Ydim_block*(this.Xdim_block - 1) - this.block_x*this.Xdim_block
				this.Neighbour_corn_id(face, B) = this.Ydim_block*(this.Xdim_block - 2) - this.block_x*this.Xdim_block
				this.Neighbour_corn_id(face, C) = this.Ydim_block*(this.Xdim_block) - this.block_x*this.Xdim_block
				this.border(face, down) = 2
				this.rot(face, down) = 2
				this.Neighbours_face(face, down) = 4
				this.Neighbours_corn_face(face, B) = 4
				this.Neighbours_corn_face(face, C) = 4
				this.Neighb_dir(face, down) = down
				this.Neighb_corn_dir(face, B) = B
				this.Neighb_corn_dir(face, C) = C

			else if ( face == 6 ) then
				this.Neighbours_face(face, down) = 2
				this.Neighbours_corn_face(face, B) = 2
				this.Neighbours_corn_face(face, C) = 2
			end if
		end if


		if (this.block_x == 0) then  ! left

			this.Neighbour_id(face, left) = this.block_y + this.Xdim_block*(this.Xdim_block - 1)
			this.Neighbour_corn_id(face, C) = this.block_y + this.Xdim_block*(this.Xdim_block - 1) - 1
			this.Neighbour_corn_id(face, D) = this.block_y + this.Xdim_block*(this.Xdim_block - 1) + 1
			this.Neighbours_face(face, left) = face - 1
			this.Neighbours_corn_face(face, C) = face - 1
			this.Neighbours_corn_face(face, D) = face - 1
			this.My_dir(face, D) = left
			this.My_dir(face, C) = left

			if ( face == 1 ) then
				this.Neighbour_id(face, left) = this.Xdim_block*this.block_y
				this.Neighbour_corn_id(face, C) = this.Xdim_block*(this.block_y - 1)
				this.Neighbour_corn_id(face, D) = this.Xdim_block*(this.block_y + 1)
				this.border(face, left) = -1
				this.rot(face, left) = -1
				this.Neighbours_face(face, left) = 5
				this.Neighbours_corn_face(face, C) = 5
				this.Neighbours_corn_face(face, D) = 5
				this.Neighb_dir(face, left) = down
				this.Neighb_corn_dir(face, D) = C
				this.Neighb_corn_dir(face, C) = B

			else if ( face == 6 ) then
				this.Neighbour_id(face, left) = this.Xdim_block*(this.Ydim_block - this.block_y) - 1
				this.Neighbour_corn_id(face, D) = this.Xdim_block*(this.Ydim_block - this.block_y - 1) - 1
				this.Neighbour_corn_id(face, C) = this.Xdim_block*(this.Ydim_block - this.block_y + 1) - 1
				this.border(face, left) = 1
				this.rot(face, left) = 1
				this.Neighbours_face(face, left) = 5
				this.Neighbours_corn_face(face, C) = 5
				this.Neighbours_corn_face(face, D) = 5
				this.Neighb_dir(face, left) = up
				this.Neighb_corn_dir(face, D) = A
				this.Neighb_corn_dir(face, C) = D

			else if ( face == 2 ) then
				this.Neighbours_face(face, left) = 5
				this.Neighbours_corn_face(face, D) = 5
				this.Neighbours_corn_face(face, C) = 5
			end if
		end if

		if (this.block_y == this.Ydim_block - 1 .and. this.block_x == this.Xdim_block - 1) this.Neighbour_corn_id (:, A) = -1
		if (this.block_y == 0 .and. this.block_x == this.Xdim_block - 1) this.Neighbour_corn_id (:, B) = -1
		if (this.block_y == 0 .and. this.block_x == 0) this.Neighbour_corn_id (:, C) = -1
		if (this.block_y == this.Ydim_block - 1 .and. this.block_x == 0) this.Neighbour_corn_id (:, D) = -1
	end do

	! print *, "My = ", id, this.Neighbour_corn_id(6, :)

End Subroutine



Subroutine Displacement(this, face, dir, displ)
	Class(parallel) :: this
	Integer(4), intent(in) :: face, dir
	Integer(4), intent(out) :: displ(1:this.step*this.Xsize)
	Integer(4) k, i, n, x

	n = this.step*this.Xsize
	x=this.Xsize + 2*this.step

	if ( dir == this.up .or. dir == this.down ) then

		select case(this.rot(face, dir) == 2)
		case (.false.)

			displ(1) = 0
			do k = 2, this.Xsize
				displ(k) = displ(k - 1) + 1
			end do

			if(this.step > 1) then
				do i = this.Xsize+1, n
					displ(i) = displ(i - this.Xsize) + x
				end do
			end if

		case (.true.)

			displ(1) = x*(this.step - 1) + this.Xsize - 1
			do k = 2, this.Xsize
				displ(k) = displ(k - 1) - 1
			end do

			if(this.step > 1) then
				do i = this.Xsize+1, n
					displ(i) = displ(i - this.Xsize) - x
				end do
			end if

		end select


	else if ( dir == this.right .or. dir == this.left ) then

		select case(this.border(face, dir))
		case (0)

				displ(1) = 0
				do k = 2, this.Ysize
					displ(k) = displ(k - 1) + x
				end do

				if(this.step > 1) then
					do i = this.Ysize+1, n
						displ(i) = displ(i - this.Ysize) + 1
					end do
				end if

		case (1)

			displ(1) = (this.Ysize - 1)*x
			do k = 2, this.Ysize
				displ(k) = displ(k - 1) - x
			end do

			if(this.step > 1) then
				do i = this.Ysize+1, n
					displ(i) = displ(i - this.Ysize) + 1
				end do
			end if

		case (-1)

			displ(1) = this.step - 1
			do k = 2, this.Ysize
				displ(k) = displ(k - 1) + x
			end do

			if(this.step > 1) then
				do i = this.Ysize+1, n
					displ(i) = displ(i - this.Ysize) - 1
				end do
			end if

		end select

	end if



End Subroutine


Subroutine Displacement_corn(this, face, dir, displ)
	Class(parallel) :: this
	Integer(4), intent(in) :: face, dir
	Integer(4), intent(out) :: displ(1:this.step*this.step)
	Integer(4) k, i, n, x

	n = this.step*this.step
	x=this.Xsize + 2*this.step

	if ( this.My_dir(face, dir) == this.right .or. this.My_dir(face, dir) == this.left ) then


		select case(this.border(face, this.My_dir(face, dir)))
		case (0)

			displ(1) = 0
			do k = 2, this.step
				displ(k) = displ(k - 1) + x
			end do

			if(this.step > 1) then
				do i = this.step+1, n
					displ(i) = displ(i - this.step) + 1
				end do
			end if

		case (1)

			displ(1) = (this.step - 1)*x
			do k = 2, this.step
				displ(k) = displ(k - 1) - x
			end do

			if(this.step > 1) then
				do i = this.step+1, n
					displ(i) = displ(i - this.step) + 1
				end do
			end if

		case (-1)

			displ(1) = this.step - 1
			do k = 2, this.step
				displ(k) = displ(k - 1) + x
			end do

			if(this.step > 1) then
				do i = this.step+1, n
					displ(i) = displ(i - this.step) - 1
				end do
			end if

		end select

	else if ( this.My_dir(face, dir) == this.up .or. this.My_dir(face, dir) == this.down ) then

		i = this.My_dir(face, dir)
		select case(this.rot(face, i) == 2)
		case (.false.)

			do k = 1, this.step
				displ(k) = k - 1
			end do

			if(this.step > 1) then
				do i = this.step+1, n
					displ(i) = displ(i - this.step) + x
				end do
			end if

		case (.true.)

			displ(1) = x*(this.step - 1) + this.step - 1
			do k = 2, this.step
				displ(k) = displ(k - 1) - 1
			end do

			if(this.step > 1) then
				do i = this.step+1, n
					displ(i) = displ(i - this.step) - x
				end do
			end if
		end select


	else
		do k = 1, this.step
			displ(k) = k - 1
		end do

		if(this.step > 1) then
			do i = this.step+1, n
				displ(i) = displ(i - this.step) + x
			end do
		end if

	end if

End Subroutine



end module