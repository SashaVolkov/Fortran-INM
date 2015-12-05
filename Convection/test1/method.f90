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
			Procedure :: to_print => met_Print
! 			Procedure :: SchemeParam => met_SchemeParam
! 			Procedure ::  BornParam=> met_func_BornParam

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
		Real(8), Intent(inout) :: trans_mass(g.ns_y: g.nf_y, g.ns_x - g.bstep : g.nf_x + g.fstep) 
		Integer(4) y

		y=g.StepsY

	if ( g.np>1 ) then

		if ( g.id==g.np-1 ) then
			call MPI_Send(trans_mass(1, g.nf_x+1-g.bstep), g.bstep*y, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, this.ier)
		else
			call MPI_Send(trans_mass(1, g.nf_x+1-g.bstep), g.bstep*y, MPI_DOUBLE_PRECISION, g.id+1, g.id+1, MPI_COMM_WORLD, this.ier)
		end if

		if ( g.id==0 ) then
			call MPI_Recv(trans_mass(1, g.ns_x - g.bstep), g.bstep*y, MPI_DOUBLE_PRECISION, g.np-1, g.id, MPI_COMM_WORLD, this.status, this.ier);
		else
			call MPI_Recv(trans_mass(1, g.ns_x - g.bstep), g.bstep*y, MPI_DOUBLE_PRECISION, g.id-1, g.id, MPI_COMM_WORLD, this.status, this.ier);
		end if


		if ( g.fstep > 0) then

			if ( g.id==0 ) then
				call MPI_Send(trans_mass(1, g.ns_x), g.fstep*y, MPI_DOUBLE_PRECISION, g.np-1, g.id, MPI_COMM_WORLD, this.ier);
			else
				call MPI_Send(trans_mass(1, g.ns_x), g.fstep*y, MPI_DOUBLE_PRECISION, g.id-1, g.id, MPI_COMM_WORLD, this.ier);
			end if

			if ( g.id==g.np-1 ) then
				call MPI_Recv(trans_mass(1, g.nf_x + 1), g.fstep*y, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, this.status, this.ier)
			else
				call MPI_Recv(trans_mass(1, g.nf_x + 1), g.fstep*y, MPI_DOUBLE_PRECISION, g.id+1, g.id+1, MPI_COMM_WORLD, this.status, this.ier)
			end if

		end if
	end if

	End Subroutine


	Subroutine met_Print(this, f, g, t, name, timeset)

		Class(met) :: this
		Class(func) :: f
		Class(grid) :: g

		Integer x, y, i, j, request(g.np)
		Integer(4), Intent(In) :: t, timeset
		character(40), Intent(In) :: name

		Real(8) W_mass(1:g.StepsY, 1:g.StepsX)

		do j=1, g.StepsX
			do i = 1, g.StepsY
				W_mass(i, j) = 0
			end do
		end do

		do j = g.ns_x, g.nf_x
			do i = 1, g.StepsY
				W_mass(i, j) = f.d(i, j)
			end do
		end do


		if ( g.np > 1 ) call MPI_Gather(W_mass(1, g.id*g.Xsize + 1), g.Xsize*g.StepsY, MPI_DOUBLE_PRECISION, W_mass, g.Xsize*g.StepsY,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, this.ier)

		if ( g.id == 0 ) then
			print *, "Writing"
			open(14,file='/home/sasha/Fortran/Convection/datFiles/'//name)
			do y = 1, g.StepsY
				write(14,'(3(e20.12))') (Real(x-1,8)*g.dx, Real(y-1,8)*g.dy, W_mass(y,x), x=1,g.StepsX)
			end do
			close(14)
		end if

	End Subroutine


! 	Subroutine met_func_BornParam(this, f, g)

! 		Class(met) :: this
! 		Class(func) :: f
! 		Class(grid) :: g

! 		Integer(4) i

! 		this.Fmax = MAXVAL(f.d(g.ns_x: g.nf_x))
! 		this.Fmin = MINVAL(f.d(g.ns_x: g.nf_x))
! 		call MPI_AllReduce(this.Fmax, this.AbsFmax_l, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, this.ier)


! 		this.Fmed = this.AbsFmax_l/2
! 		this.w_counter = 0

! 		if ( this.Fmax < this.Fmed ) then
! 			!nothing
! 		elseif ( this.Fmin > this.Fmed ) then
! 			!all
! 			this.w_counter = (g.nf_x - g.ns_x)*g.dx
! 		elseif ( this.Fmax > this.Fmed .AND. this.Fmin < this.Fmed ) then
! 			!cycle
! 			do i = g.ns_x, g.nf_x
! 				if ( f.d(i) > this.Fmed ) then
! 					this.w_counter = this.w_counter + g.dx
! 				end if
! 			end do
! 		end if

! 		call MPI_Reduce(this.w_counter, this.width_last, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, this.ier)


! 		this.N1_mass = 0
! 		do i = g.ns_x, g.nf_x
! 			this.N1_mass = this.N1_mass + g.dx*((f.d(i-1)+f.d(i))/2.0)
! 		end do
! 		call MPI_ALLREDUCE(this.N1_mass, this.N1_All_mass_l, g.np, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, this.ier)


! 		this.N2_mass = 0
! 		do i = g.ns_x, g.nf_x
! 			this.N2_mass = this.N2_mass + g.dx*(((f.d(i-1)+f.d(i))*(f.d(i-1)+f.d(i)))/4.0)
! 		end do
! 		call MPI_ALLREDUCE(this.N2_mass, this.N2_All_mass_l, g.np, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, this.ier)


! 	End Subroutine



! 	Subroutine met_SchemeParam(this, f, g, name)

! 		Class(met) :: this
! 		Class(func) :: f
! 		Class(grid) :: g

! 		character(20), Intent(In) :: name

! 		Real(8) :: width_first, AbsFmax_f, N1_All_mass_f, N2_All_mass_f, N1_Mass_Def, N2_Mass_Def
! 		Real(8) :: AbsF, PulseStretching


! 		width_first = this.width_last
! 		AbsFmax_f = this.AbsFmax_l
! 		N1_All_mass_f = this.N1_All_mass_l
! 		N2_All_mass_f = this.N2_All_mass_l

! 		call MPI_Barrier(MPI_COMM_WORLD, this.ier)
! 		if ( g.id == 0 ) then
! 			print *, "SchemeParam start"
! 		end if

! 		call this.BornParam(f, g)
! 		call MPI_Barrier(MPI_COMM_WORLD, this.ier)

! 		if ( g.id == 0 ) then

! 			PulseStretching = width_first - this.width_last
! 			AbsF = AbsFmax_f - this.AbsFmax_l
! 			N1_Mass_Def = N1_All_mass_f - this.N1_All_mass_l
! 			N2_Mass_Def = N2_All_mass_f - this.N2_All_mass_l

! ! 			open(10,file='SchemeParam'//name, position="append")
! ! 			write(10,*) (name)
! ! 			write(10,*) ("Number of steps", g.StepsX)
! ! 			write(10,*) ("Decreasing = ", AbsF)
! ! 			write(10,*) ("Stretching = ", PulseStretching)
! ! 			write(10,*) ("N1 Mass Defect = ", N1_Mass_Def)
! ! 			write(10,*) ("N2 Mass Defect = ", N2_Mass_Def)


! 			if ( this.newfile == 1 ) then
! 				open(10,file='/home/sasha/Fortran/Convection/MethodParams/SchemeParam'//name)
! 				write(10,*) ('Steps; Decreasing; Stretching; L1; L2') !,'Decreasing','Stretching','L1','L2')
! 				write(10,FMT="(I14, A, E14.7, A, E14.7, A, E14.7, A, E14.7)") g.StepsX,";",AbsF,";", PulseStretching,";", N1_Mass_Def,";", N2_Mass_Def
! 				close(10)

! 			else
! 				open(10,file='/home/sasha/Fortran/Convection/MethodParams/SchemeParam'//name, position="append")
! ! 				write(10,*) ('Steps; Decreasing; Stretching; L1; L2') !,'Decreasing','Stretching','L1','L2')
! 				write(10,FMT="(I14, A, E14.7, A, E14.7, A, E14.7, A, E14.7)") g.StepsX,";",AbsF,";", PulseStretching,";", N1_Mass_Def,";", N2_Mass_Def
! 				close(10)

! 			end if

! 		end if


! 	End Subroutine


! 	Subroutine met_deinit(this)
! 		Class(met) :: this

! 		if (Allocated(this.k1)) Deallocate(this.k1)
! 		if (Allocated(this.k2)) Deallocate(this.k2)
! 		if (Allocated(this.k3)) Deallocate(this.k3)
! 		if (Allocated(this.k4)) Deallocate(this.k4)

! 	End Subroutine


End Module