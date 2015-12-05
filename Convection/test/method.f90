Module method

	Use modnet
	Use modfunc
	Use MPI

IMPLICIT NONE

	Public :: met

	Type met

		Integer(4) ier, newfile
		Integer(4) :: status(MPI_STATUS_SIZE)
		Real(8) Fmax, Fmin, Fmed, AbsF, AbsFmax_l, PulseStretching, width_last, param
		Real(8) w_counter, temp, N1_mass, N1_All_mass_l, N2_mass, N2_All_mass_l

		CONTAINS

			Procedure :: Message => met_Message

			Procedure :: init => met_init
! 			Procedure :: deinit => met_deinit

	End Type

	CONTAINS

	Subroutine met_init(this, g, status, ier, newfile)

		Class(met) :: this
		Class(grid) :: g
		
		Integer(4), Intent(In) :: ier, newfile
		Integer(4), Intent(In) :: status(MPI_STATUS_SIZE)

		this.ier=ier
		this.status=status
		this.newfile=newfile

	End Subroutine


	Subroutine met_Message(this, trans_mass, g)

		Class(met) :: this
		Class(grid) :: g
		Real(8), Intent(inout) :: trans_mass(g.ns - g.bstep : g.nf + g.fstep)

	if ( g.np>1 ) then

		if ( g.id==g.np-1 ) then
			call MPI_Send(trans_mass(g.nf+1-g.bstep), g.bstep, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, this.ier)
		else
			call MPI_Send(trans_mass(g.nf+1-g.bstep), g.bstep, MPI_DOUBLE_PRECISION, g.id+1, g.id+1, MPI_COMM_WORLD, this.ier)
		end if

		if ( g.id==0 ) then
			call MPI_Recv(trans_mass(g.ns - g.bstep), g.bstep, MPI_DOUBLE_PRECISION, g.np-1, g.id, MPI_COMM_WORLD, this.status, this.ier);
		else
			call MPI_Recv(trans_mass(g.ns - g.bstep), g.bstep, MPI_DOUBLE_PRECISION, g.id-1, g.id, MPI_COMM_WORLD, this.status, this.ier);
		end if


		if ( g.fstep > 0) then

			if ( g.id==0 ) then
				call MPI_Send(trans_mass(g.ns), g.fstep, MPI_DOUBLE_PRECISION, g.np-1, g.id, MPI_COMM_WORLD, this.ier);
			else
				call MPI_Send(trans_mass(g.ns), g.fstep, MPI_DOUBLE_PRECISION, g.id-1, g.id, MPI_COMM_WORLD, this.ier);
			end if

			if ( g.id==g.np-1 ) then
				call MPI_Recv(trans_mass(g.nf + 1), g.fstep, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, this.status, this.ier)
			else
				call MPI_Recv(trans_mass(g.nf + 1), g.fstep, MPI_DOUBLE_PRECISION, g.id+1, g.id+1, MPI_COMM_WORLD, this.status, this.ier)
			end if

		end if
	end if

	End Subroutine

End Module