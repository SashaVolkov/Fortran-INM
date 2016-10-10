module messenger

	use parallel_cubic, Only: parallel
	use func_var, Only: f_var
	use mpi

	implicit none


	Private
	Public :: message

	Type message

	integer(4) :: var_count = 3
	integer(4) snd_stat(MPI_STATUS_SIZE, 3, 4, 6), rcv_stat(MPI_STATUS_SIZE, 3, 4, 6)
	integer(4) ier, np, snd_req(3, 4, 6), rcv_req(3, 4, 6), comm(3), grid_type, vec_only

		CONTAINS
		Procedure, Public :: msg => msg
		Procedure, Public :: init => init
		Procedure, Private :: Simple_msg => Simple_msg
		Procedure, Private :: Waiter => Waiter
		Procedure, Private :: Halo_Rotations => Halo_Rotations
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



subroutine msg(this, f, paral, vec_only)
	Class(message) :: this
	Class(f_var) :: f
	Class(parallel) :: paral
	integer(4), intent(in) :: vec_only

	this.vec_only = vec_only

	call this.Simple_msg(paral, f)

	call this.Waiter(paral)

	if (this.grid_type == 0) then
		call this.Halo_Rotations(paral, f.lon_vel, 1, 0)
		call this.Halo_Rotations(paral, f.lat_vel, 0, 1)
	end if


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

			if (this.vec_only == 0) then
				call MPI_IRecv(f.h_height(rx, ry, face), 1, paral.halo(face, i), neib_id, rcv_tag, this.comm(1), this.rcv_req(1, i, face), this.ier)
				call MPI_ISend(f.h_height(sx, sy, face), 1, paral.halo(face, i), neib_id, snd_tag, this.comm(1), this.snd_req(1, i, face), this.ier)
			else
				if(paral.border(face, i) == 0 .or. paral.border(face, i) == 2  .or. this.grid_type == 1) then ! coordinate rotation from one face to other if rotation 90 or -90 deg
					call MPI_IRecv(f.lon_vel(rx, ry, face), 1, paral.halo(face, i), neib_id, rcv_tag, this.comm(2), this.rcv_req(2, i, face), this.ier) ! x -> x
					call MPI_IRecv(f.lat_vel(rx, ry, face), 1, paral.halo(face, i), neib_id, rcv_tag, this.comm(3), this.rcv_req(3, i, face), this.ier) ! y -> y
				else
					call MPI_IRecv(f.lon_vel(rx, ry, face), 1, paral.halo(face, i), neib_id, rcv_tag, this.comm(3), this.rcv_req(3, i, face), this.ier) ! x -> y
					call MPI_IRecv(f.lat_vel(rx, ry, face), 1, paral.halo(face, i), neib_id, rcv_tag, this.comm(2), this.rcv_req(2, i, face), this.ier) ! y -> x
				end if

				call MPI_ISend(f.lon_vel(sx, sy, face), 1, paral.halo(face, i), neib_id, snd_tag, this.comm(2), this.snd_req(2, i, face), this.ier)
				call MPI_ISend(f.lat_vel(sx, sy, face), 1, paral.halo(face, i), neib_id, snd_tag, this.comm(3), this.snd_req(3, i, face), this.ier)
			end if
		end do
	end do

end subroutine



subroutine Waiter(this, paral)

	Class(message) :: this
	Class(parallel) :: paral

	if (this.vec_only == 0) then
		call MPI_Waitall(6*4, this.rcv_req(1, :, :), this.rcv_stat, this.ier)
		call MPI_Waitall(6*4, this.snd_req(1, :, :), this.snd_stat, this.ier)
	else
		call MPI_Waitall(6*4*2, this.rcv_req(2:3, :, :), this.rcv_stat, this.ier)
		call MPI_Waitall(6*4*2, this.snd_req(2:3, :, :), this.snd_stat, this.ier)
	end if

end subroutine


subroutine Halo_Rotations(this, paral, func_mass, x_vec, y_vec)

	Class(message) :: this
	Class(parallel) :: paral

	real(8), intent(inout) :: func_mass(paral.first_x:paral.last_x, paral.first_y:paral.last_y, 6)
	integer(4), intent(in) :: x_vec, y_vec
	integer(4) face, up, right, left, down, i, ns_xy(2), nf_xy(2), k

	up = 1; right = 2; down = 3; left = 4

	do face = 1, 6
		do i = 1, 4

			ns_xy(:) = paral.rcv_xy(face, i, :)
			if(i == up .or. i == down) then
				nf_xy(1) = ns_xy(1) + paral.Xsize - 1
				nf_xy(2) = ns_xy(2) + paral.step - 1
			else
				nf_xy(1) = ns_xy(1) + paral.step - 1
				nf_xy(2) = ns_xy(2) + paral.Ysize - 1
			end if

			select case(paral.border(face, i) == 2)
				case(.false.)
					if(y_vec == 1 .and. paral.rot(face, i) == 1) func_mass(ns_xy(1):nf_xy(1), ns_xy(2):nf_xy(2), face) = - func_mass(ns_xy(1):nf_xy(1), ns_xy(2):nf_xy(2), face)
					if(x_vec == 1 .and. paral.rot(face, i) == -1) func_mass(ns_xy(1):nf_xy(1), ns_xy(2):nf_xy(2), face) = - func_mass(ns_xy(1):nf_xy(1), ns_xy(2):nf_xy(2), face)
				case(.true.)
					if(x_vec == 1 .or. y_vec == 1) func_mass(ns_xy(1):nf_xy(1), ns_xy(2):nf_xy(2), face) = - func_mass(ns_xy(1):nf_xy(1), ns_xy(2):nf_xy(2), face)
			end select

		end do
	end do

end subroutine



end module