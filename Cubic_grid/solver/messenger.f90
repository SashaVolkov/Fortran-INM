module messenger

	use parallel_cubic, Only: parallel
	use func_var, Only: f_var
	use mpi

	implicit none


	Private
	Public :: message

	Type message

	integer(4) :: var_count = 3
	integer(4) snd_stat(MPI_STATUS_SIZE, 4, 6, 3), rcv_stat(MPI_STATUS_SIZE, 4, 6, 3)
	integer(4) ier, np, snd_req(4, 6, 3), rcv_req(4, 6, 3), comm(3), grid_type, vec_only

		CONTAINS
		Procedure, Public :: msg => msg
		Procedure, Public :: init => init
		Procedure, Private :: Simple_msg => Simple_msg
		Procedure, Private :: Waiter => Waiter
	End Type


	CONTAINS



subroutine init(this, grid_type)

	Class(message) :: this
	integer(4), intent(in) :: grid_type
	integer(4) i

	this.grid_type = grid_type

	call MPI_Comm_size(MPI_COMM_WORLD,this.np,this.ier)

	do i = 1, this.var_count
		call MPI_COMM_DUP(MPI_COMM_WORLD, this.comm(i), this.ier)
	end do

end subroutine



subroutine msg(this, f, paral)
	Class(message) :: this
	Class(f_var) :: f
	Class(parallel) :: paral

	call this.Simple_msg(paral, f)
	call this.Waiter(paral)


end subroutine




subroutine Simple_msg(this, paral, f)

	Class(message) :: this
	Class(f_var) :: f
	Class(parallel) :: paral
	integer(4) i, face, rcv_tag, snd_tag
	integer(4) rx, ry, sx, sy, neib_id


	do face = 1, 6
		do i = 1, 4

			neib_id = paral.Neighbour_id(face, i)
			rx = paral.rcv_xy(face, i, 1);  ry = paral.rcv_xy(face, i, 2)
			sx = paral.snd_xy(face, i, 1);  sy = paral.snd_xy(face, i, 2)
			snd_tag = (neib_id + 1)*6*4 + paral.Neighbours_face(face, i)*4 + paral.Neighb_dir(face, i)
			rcv_tag = (paral.id + 1)*6*4  + face*4 + i

				call MPI_IRecv(f.h_height(rx, ry, face), 1, paral.halo(face, i), neib_id, rcv_tag, this.comm(1), this.rcv_req(i, face, 1), this.ier)
				call MPI_IRecv(f.lon_vel(rx, ry, face), 1, paral.halo(face, i), neib_id, rcv_tag, this.comm(2), this.rcv_req(i, face, 2), this.ier) ! x -> x
				call MPI_IRecv(f.lat_vel(rx, ry, face), 1, paral.halo(face, i), neib_id, rcv_tag, this.comm(3), this.rcv_req(i, face, 3), this.ier) ! y -> y

				call MPI_ISend(f.h_height(sx, sy, face), 1, paral.halo(face, i), neib_id, snd_tag, this.comm(1), this.snd_req(i, face, 1), this.ier)
				call MPI_ISend(f.lon_vel(sx, sy, face), 1, paral.halo(face, i), neib_id, snd_tag, this.comm(2), this.snd_req(i, face, 2), this.ier)
				call MPI_ISend(f.lat_vel(sx, sy, face), 1, paral.halo(face, i), neib_id, snd_tag, this.comm(3), this.snd_req(i, face, 3), this.ier)
		end do
	end do

end subroutine



subroutine Waiter(this, paral)

	Class(message) :: this
	Class(parallel) :: paral

		call MPI_Waitall(6*4*3, this.rcv_req(:, :, :), this.rcv_stat, this.ier)
		call MPI_Waitall(6*4*3, this.snd_req(:, :, :), this.snd_stat, this.ier)

end subroutine



end module