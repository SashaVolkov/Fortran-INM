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
		! Procedure, Private :: Rot_minus_90_msg => Rot_minus_90_msg
		! Procedure, Private :: Rot_90_msg => Rot_90_msg
		! Procedure, Private :: Rot_180_msg => Rot_180_msg
		Procedure, Private :: Waiter => Waiter
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

	integer(4) face, i


	call this.Simple_msg(paral, f)

	call this.Waiter(paral)

end subroutine




subroutine Simple_msg(this, paral, f)

	Class(message) :: this
	Class(f_var) :: f
	Class(parallel) :: paral
	integer(4) i, face
	integer(4) rx, ry, sx, sy, neib_id

	do face = 1, 6
		do i = 1, 4


	neib_id = paral.Neighbour_id(face, i)
	rx = paral.rcv_xy(i, 1);  ry = paral.rcv_xy(i, 2)
	sx = paral.snd_xy(i, 1);  sy = paral.snd_xy(i, 2)


			call MPI_IRecv(f.h_height(rx, ry, face), 1, paral.grey(i), neib_id, paral.id, this.comm(1), this.rcv_req(1, i, face), this.ier)
			call MPI_IRecv(f.x_vel(rx, ry, face), 1, paral.grey(i), neib_id, paral.id, this.comm(2), this.rcv_req(2, i, face), this.ier)
			call MPI_IRecv(f.y_vel(rx, ry, face), 1, paral.grey(i), neib_id, paral.id, this.comm(3), this.rcv_req(3, i, face), this.ier)

			call MPI_ISend(f.h_height(sx, sy, face), 1, paral.grey(i), neib_id, neib_id, this.comm(1), this.snd_req(1, i, face), this.ier)
			call MPI_ISend(f.x_vel(sx, sy, face), 1, paral.grey(i), neib_id, neib_id, this.comm(2), this.snd_req(2, i, face), this.ier)
			call MPI_ISend(f.y_vel(sx, sy, face), 1, paral.grey(i), neib_id, neib_id, this.comm(3), this.snd_req(3, i, face), this.ier)

		end do
	end do

end subroutine



subroutine Waiter(this, paral)

	Class(message) :: this
	Class(parallel) :: paral
	integer(4) i, face, FLAG


	call MPI_Waitall(6*4*3, this.rcv_req, this.rcv_stat, this.ier)
	call MPI_Waitall(6*4*3, this.snd_req, this.snd_stat, this.ier)


end subroutine




end module