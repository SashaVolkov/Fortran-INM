module messenger

	use parallel_cubic, Only: parallel
	use mpi

	implicit none


	Private
	Public :: message

	Type message

	integer(4) :: var_count = 3
	integer(4) snd_stat(MPI_STATUS_SIZE, 4, 6, 3), rcv_stat(MPI_STATUS_SIZE, 4, 6, 3)
	integer(4) ier, np, snd_req(4, 6, 3), rcv_req(4, 6, 3), comm(3), grid_type, vec_only
	integer(4) snd_xy(6, 4, 2), rcv_xy(6, 4, 2), halo(6, 4), first_x, first_y, last_x, last_y
	integer(4) Neighbour_id(1:6, 1:4), border(6, 4), Neighbours_face(6, 4), id, Neighb_dir(6,4)

		CONTAINS
		Procedure, Public :: msg => msg
		Procedure, Public :: init => init
		Procedure, Private :: Simple_msg => Simple_msg
		Procedure, Private :: Waiter => Waiter
	End Type


	CONTAINS



subroutine init(this, grid_type, paral)

	Class(message) :: this
	Class(parallel) :: paral
	integer(4), intent(in) :: grid_type
	integer(4) i

	this.grid_type = grid_type;  this.id = paral.id
	this.halo(:,:) = paral.halo(:,:);  this.Neighbours_face(:,:) = paral.Neighbours_face(:,:)
	this.Neighb_dir(:,:) = paral.Neighb_dir(:,:);  this.Neighbour_id(:,:) = paral.Neighbour_id(:,:)
	this.rcv_xy(:,:,:) = paral.rcv_xy(:,:,:);  this.snd_xy(:,:,:) = paral.snd_xy(:,:,:)

	this.first_x = paral.first_x;  this.first_y = paral.first_y
	this.last_x = paral.last_x;  this.last_y = paral.last_y

	call MPI_Comm_size(MPI_COMM_WORLD,this.np,this.ier)

	do i = 1, this.var_count
		call MPI_COMM_DUP(MPI_COMM_WORLD, this.comm(i), this.ier)
	end do

end subroutine



subroutine msg(this, level, lon_vel, lat_vel)
	Class(message) :: this
	Real(8), Intent(in) :: level(this.first_x:this.last_x, this.first_y:this.last_y, 6), lon_vel(this.first_x:this.last_x, this.first_y:this.last_y, 6), lat_vel(this.first_x:this.last_x, this.first_y:this.last_y, 6)

	call this.Simple_msg(level, lon_vel, lat_vel)
	call this.Waiter()


end subroutine




subroutine Simple_msg(this, level, lon_vel, lat_vel)

	Class(message) :: this
	Real(8), Intent(in) :: level(this.first_x:this.last_x, this.first_y:this.last_y, 6), lon_vel(this.first_x:this.last_x, this.first_y:this.last_y, 6), lat_vel(this.first_x:this.last_x, this.first_y:this.last_y, 6)
	integer(4) i, face, rcv_tag, snd_tag
	integer(4) rx, ry, sx, sy, neib_id


	do face = 1, 6
		do i = 1, 4

			neib_id = this.Neighbour_id(face, i)
			rx = this.rcv_xy(face, i, 1);  ry = this.rcv_xy(face, i, 2)
			sx = this.snd_xy(face, i, 1);  sy = this.snd_xy(face, i, 2)
			snd_tag = (neib_id + 1)*6*4 + this.Neighbours_face(face, i)*4 + this.Neighb_dir(face, i)
			rcv_tag = (this.id + 1)*6*4  + face*4 + i

				call MPI_IRecv(level(rx, ry, face), 1, this.halo(face, i), neib_id, rcv_tag, this.comm(1), this.rcv_req(i, face, 1), this.ier)
				call MPI_IRecv(lon_vel(rx, ry, face), 1, this.halo(face, i), neib_id, rcv_tag, this.comm(2), this.rcv_req(i, face, 2), this.ier) ! x -> x
				call MPI_IRecv(lat_vel(rx, ry, face), 1, this.halo(face, i), neib_id, rcv_tag, this.comm(3), this.rcv_req(i, face, 3), this.ier) ! y -> y

				call MPI_ISend(level(sx, sy, face), 1, this.halo(face, i), neib_id, snd_tag, this.comm(1), this.snd_req(i, face, 1), this.ier)
				call MPI_ISend(lon_vel(sx, sy, face), 1, this.halo(face, i), neib_id, snd_tag, this.comm(2), this.snd_req(i, face, 2), this.ier)
				call MPI_ISend(lat_vel(sx, sy, face), 1, this.halo(face, i), neib_id, snd_tag, this.comm(3), this.snd_req(i, face, 3), this.ier)
		end do
	end do

end subroutine



subroutine Waiter(this)

	Class(message) :: this

		call MPI_Waitall(6*4*3, this.rcv_req(:, :, :), this.rcv_stat, this.ier)
		call MPI_Waitall(6*4*3, this.snd_req(:, :, :), this.snd_stat, this.ier)

end subroutine



end module