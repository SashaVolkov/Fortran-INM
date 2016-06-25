module messenger

	use parallel_cubic, Only: parallel
	use func_var, Only: f_var

	implicit none

	include"mpif.h"

	Private
	Public :: message

	Type message

	integer(4) :: var_count = 3
	integer(4) snd_stat(MPI_STATUS_SIZE, 3, 4, 6), rcv_stat(MPI_STATUS_SIZE, 3, 4, 6)
	integer(4) ier, snd_req(3, 4, 6), rcv_req(3, 4, 6), comm(3)

		CONTAINS
		Procedure, Public :: msg => msg
		Procedure, Public :: init => init
		Procedure, Private :: Simple_msg => Simple_msg
		Procedure, Private :: Waiter => Waiter
		Procedure, Private :: Halo_Rotations => Halo_Rotations
		Procedure, Private :: swap => swap
	End Type


	CONTAINS



subroutine init(this)

	Class(message) :: this
	integer(4) i

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

! 	call this.Halo_Rotations(paral, f.h_height, 0, 0)
	! call this.Halo_Rotations(paral, f.x_vel, 1, 0)
	! call this.Halo_Rotations(paral, f.y_vel, 0, 1)

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
			rcv_tag = (paral.id + 1) * 6 + face
			snd_tag = (neib_id + 1) * 6 + paral.Neighbours_face(face, i)

			call MPI_IRecv(f.h_height(rx, ry, face), 1, paral.halo(face, i), neib_id, rcv_tag, this.comm(1), this.rcv_req(1, i, face), this.ier)

			if(paral.border(face, i) == 0 .or. paral.border(face, i) == 2) then ! coordinate rotation from one face to other if rotation 90 or -90 deg
				call MPI_IRecv(f.x_vel(rx, ry, face), 1, paral.halo(face, i), neib_id, rcv_tag, this.comm(2), this.rcv_req(2, i, face), this.ier) ! x -> x
				call MPI_IRecv(f.y_vel(rx, ry, face), 1, paral.halo(face, i), neib_id, rcv_tag, this.comm(3), this.rcv_req(3, i, face), this.ier) ! y -> y
			else
				call MPI_IRecv(f.x_vel(rx, ry, face), 1, paral.halo(face, i), neib_id, rcv_tag, this.comm(3), this.rcv_req(3, i, face), this.ier) ! x -> y
				call MPI_IRecv(f.y_vel(rx, ry, face), 1, paral.halo(face, i), neib_id, rcv_tag, this.comm(2), this.rcv_req(2, i, face), this.ier) ! y -> x
			end if

			call MPI_ISend(f.h_height(sx, sy, face), 1, paral.halo(face, i), neib_id, snd_tag, this.comm(1), this.snd_req(1, i, face), this.ier)
			call MPI_ISend(f.x_vel(sx, sy, face), 1, paral.halo(face, i), neib_id, snd_tag, this.comm(2), this.snd_req(2, i, face), this.ier)
			call MPI_ISend(f.y_vel(sx, sy, face), 1, paral.halo(face, i), neib_id, snd_tag, this.comm(3), this.snd_req(3, i, face), this.ier)

		end do
	end do

end subroutine



subroutine Waiter(this, paral)

	Class(message) :: this
	Class(parallel) :: paral

	call MPI_Waitall(6*4*this.var_count, this.rcv_req, this.rcv_stat, this.ier)
	call MPI_Waitall(6*4*this.var_count, this.snd_req, this.snd_stat, this.ier)

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

			select case(paral.border(face, i))
				case(-1) ! step rot


					if(i == up .or. i == down) then
						do k = ns_xy(1), nf_xy(1), 1
							call this.swap(func_mass(k, ns_xy(2):nf_xy(2), face), paral.step)
						end do
					else
						do k = ns_xy(2), nf_xy(2), 1
							call this.swap(func_mass(ns_xy(1):nf_xy(1), k, face), paral.step)
						end do
					end if


					if(y_vec == 1 .and. paral.rot(face, i) == 1) func_mass(ns_xy(1):nf_xy(1), ns_xy(2):nf_xy(2), face) = - func_mass(ns_xy(1):nf_xy(1), ns_xy(2):nf_xy(2), face)
					if(x_vec == 1 .and. paral.rot(face, i) == -1) func_mass(ns_xy(1):nf_xy(1), ns_xy(2):nf_xy(2), face) = - func_mass(ns_xy(1):nf_xy(1), ns_xy(2):nf_xy(2), face)


				case(1) ! side rot

				if(i == up .or. i == down) then
					do k = ns_xy(2), nf_xy(2), 1
						call this.swap(func_mass(ns_xy(1):nf_xy(1), k, face), paral.Xsize)
					end do
				else
					do k = ns_xy(1), nf_xy(1), 1
						call this.swap(func_mass(k, ns_xy(2):nf_xy(2), face), paral.Ysize)
					end do
				end if

				if(y_vec == 1 .and. paral.rot(face, i) == 1) func_mass(ns_xy(1):nf_xy(1), ns_xy(2):nf_xy(2), face) = - func_mass(ns_xy(1):nf_xy(1), ns_xy(2):nf_xy(2), face)
				if(x_vec == 1 .and. paral.rot(face, i) == -1) func_mass(ns_xy(1):nf_xy(1), ns_xy(2):nf_xy(2), face) = - func_mass(ns_xy(1):nf_xy(1), ns_xy(2):nf_xy(2), face)

					! print *, face, i, paral.id

				case(2) ! full rot

				if(i == up .or. i == down) then
					do k = ns_xy(2), nf_xy(2), 1
						call this.swap(func_mass(ns_xy(1):nf_xy(1), k, face), paral.Xsize)
					end do
					do k = ns_xy(1), nf_xy(1), 1
						call this.swap(func_mass(k, ns_xy(2):nf_xy(2), face), paral.step)
					end do
				else
					do k = ns_xy(1), nf_xy(1), 1
						call this.swap(func_mass(k, ns_xy(2):nf_xy(2), face), paral.Ysize)
					end do
					do k = ns_xy(2), nf_xy(2), 1
						call this.swap(func_mass(ns_xy(1):nf_xy(1), k, face), paral.step)
					end do
				end if
				if(x_vec ==1 .or. y_vec == 1) func_mass(ns_xy(1):nf_xy(1), ns_xy(2):nf_xy(2), face) = - func_mass(ns_xy(1):nf_xy(1), ns_xy(2):nf_xy(2), face)

				case default

			end select



		end do
	end do

end subroutine



subroutine swap(this, mass, n)

	Class(message) :: this
	integer(4), intent(in) :: n
	real(8), intent(inout) :: mass(1:n)
	real(8) temp
	integer(4) i
	i = 0

	do while(1+i < n-i)
		temp = mass(1+i)
		mass(1+i) = mass(n-i)
		mass(n-i) = temp

		i = i + 1
	end do

end subroutine



end module