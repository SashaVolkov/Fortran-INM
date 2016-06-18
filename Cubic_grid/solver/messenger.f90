module messenger

	use parallel_cubic, Only: parallel
	use func_var, Only: f_var

	implicit none

	include"mpif.h"

	Private
	Public :: message

	Type message

	integer(4) stat(MPI_STATUS_SIZE, 6, 4, 3), snd_req(6, 4, 3), rcv_req(6, 4, 3)
	integer(4) ier

		CONTAINS
		Procedure, Public :: msg=> msg
		Procedure, Private :: Simple_msg => Simple_msg
		Procedure, Private :: Rot_minus_90_msg => Rot_minus_90_msg
		Procedure, Private :: Rot_90_msg => Rot_90_msg
		Procedure, Private :: Rot_180_msg => Rot_180_msg
	End Type


	CONTAINS



subroutine msg(this, f, paral)
	Class(message) :: this
	Class(f_var) :: f
	Class(parallel) :: paral

	integer(4) face, i

	do face = 1, 1
		do i = 1, 4

	select case(paral.border(face, i))
		case(0)
			call this.Simple_msg(paral, f, face, i)
			print *, " 0", i, paral.Neighbours_face(face, i), paral.id, paral.Neighbour_id(face, i)
		case(-1)
			call this.Rot_minus_90_msg()
			print *, "-1", i, paral.Neighbours_face(face, i), paral.id, paral.Neighbour_id(face, i)
		case(1)
			call this.Rot_90_msg()
			print *, " 1", i, paral.Neighbours_face(face, i), paral.id, paral.Neighbour_id(face, i)
		case(2)
			call this.Rot_180_msg()
			print *, " 2", i, paral.Neighbours_face(face, i), paral.id, paral.Neighbour_id(face, i)
		case default
			print *, "Error", paral.border(face, i), paral.id, face
	end select

		end do
	end do

end subroutine




subroutine Simple_msg(this, paral, f, face, i)

	Class(message) :: this
	Class(f_var) :: f
	Class(parallel) :: paral
	integer(4), intent(in) :: i, face
	integer(4) mp_cw, rx, ry, sx, sy, neib_id

	mp_cw = MPI_COMM_WORLD; neib_id = paral.Neighbour_id(face, i)
	rx = paral.rcv_xy(i, 1);  ry = paral.rcv_xy(i, 2)
	sx = paral.snd_xy(i, 1);  sy = paral.snd_xy(i, 2)

	call MPI_IRecv(f.h_height(face, rx, ry), 1, paral.grey(i), neib_id, paral.id, mp_cw, this.rcv_req(face, i, 1), this.ier)
	call MPI_IRecv(f.u_vel(face, rx, ry), 1, paral.grey(i), neib_id, paral.id, mp_cw, this.rcv_req(face, i, 2), this.ier)
	call MPI_IRecv(f.v_vel(face, rx, ry), 1, paral.grey(i), neib_id, paral.id, mp_cw, this.rcv_req(face, i, 3), this.ier)

	call MPI_ISend(f.h_height(face, sx, sy), 1, paral.grey(i), paral.id, neib_id, mp_cw, this.snd_req(face, i, 1), this.ier)
	call MPI_ISend(f.u_vel(face, sx, sy), 1, paral.grey(i), paral.id, neib_id, mp_cw, this.snd_req(face, i, 2), this.ier)
	call MPI_ISend(f.v_vel(face, sx, sy), 1, paral.grey(i), paral.id, neib_id, mp_cw, this.snd_req(face, i, 3), this.ier)

end subroutine




subroutine Rot_minus_90_msg(this)

	Class(message) :: this

end subroutine




subroutine Rot_90_msg(this)

	Class(message) :: this

end subroutine




subroutine Rot_180_msg(this)

	Class(message) :: this

end subroutine




end module