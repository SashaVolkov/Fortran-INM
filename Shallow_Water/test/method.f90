Module method

	Use modnet
	Use modfunc

! 	Use netcdf

IMPLICIT NONE

	Private
	Public :: met
	include"mpif.h"

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


	Subroutine met_Message(this, f, g)

		Class(met) :: this
		Class(grid) :: g
		Class(func) :: f
! 		Real(8), Intent(inout) :: trans_mass(g.first_y : g.last_y, g.first_x : g.last_x) 
		Integer(4) y, x, rc, u, d, r, l, i, j, mp_cw, mp_dp, b_grey_zone, f_grey_zone
		Integer(4) reqsl(3), reqsr(3), reqsu(3), reqsd(3), numb_oper
		integer(4) statsl(MPI_STATUS_SIZE, 3), statsr(MPI_STATUS_SIZE, 3), statsu(MPI_STATUS_SIZE, 3), statsd(MPI_STATUS_SIZE, 3)
		integer(4) s_reqsl(3), s_reqsr(3), s_reqsu(3), s_reqsd(3)
		integer(4) s_statsl(MPI_STATUS_SIZE, 3), s_statsr(MPI_STATUS_SIZE, 3), s_statsu(MPI_STATUS_SIZE, 3), s_statsd(MPI_STATUS_SIZE, 3)

		u = g.Neighb_up; d = g.Neighb_down
		r = g.Neighb_right; l = g.Neighb_left
		mp_cw = MPI_COMM_WORLD; mp_dp = MPI_DOUBLE_PRECISION
		numb_oper = 0

		y=g.Ysize + g.bstep + g.fstep
		x=g.Xsize + g.bstep + g.fstep

		call MPI_TYPE_VECTOR(x, g.bstep, y, mp_dp, b_grey_zone, this.ier)
		call MPI_TYPE_VECTOR(x, g.fstep, y, mp_dp, f_grey_zone, this.ier)

		call MPI_TYPE_COMMIT(b_grey_zone, this.ier)
		call MPI_TYPE_COMMIT(f_grey_zone, this.ier)

if (g.np > 1) then

	if (l > -1) then
		call MPI_IRecv(f.d(g.first_y, g.first_x), g.bstep*y, mp_dp, l, g.id, mp_cw, reqsl(1), this.ier);
		call MPI_IRecv(f.du(g.first_y, g.first_x), g.bstep*y, mp_dp, l, g.id, mp_cw, reqsl(2), this.ier);
		call MPI_IRecv(f.dv(g.first_y, g.first_x), g.bstep*y, mp_dp, l, g.id, mp_cw, reqsl(3), this.ier);
	end if
	if ( r > -1) then
		call MPI_IRecv(f.d(g.first_y, g.nf_x + 1), g.fstep*y, mp_dp, r, r, mp_cw, reqsr(1), this.ier);
		call MPI_IRecv(f.du(g.first_y, g.nf_x + 1), g.fstep*y, mp_dp, r, r, mp_cw, reqsr(2), this.ier);
		call MPI_IRecv(f.dv(g.first_y, g.nf_x + 1), g.fstep*y, mp_dp, r, r, mp_cw, reqsr(3), this.ier);
	end if

	if ( u > -1 ) then
		call MPI_IRecv(f.d(g.first_y, g.first_x), 1, b_grey_zone, u, g.id, mp_cw, reqsu(1), this.ier)
		call MPI_IRecv(f.du(g.first_y, g.first_x), 1, b_grey_zone, u, g.id, mp_cw, reqsu(2), this.ier)
		call MPI_IRecv(f.dv(g.first_y, g.first_x), 1, b_grey_zone, u, g.id, mp_cw, reqsu(3), this.ier)
	end if
	if ( d > -1 ) then
		call MPI_IRecv(f.d(g.nf_y+1, g.first_x), 1, f_grey_zone, d, g.id, mp_cw, reqsd(1), this.ier)
		call MPI_IRecv(f.du(g.nf_y+1, g.first_x), 1, f_grey_zone, d, g.id, mp_cw, reqsd(2), this.ier)
		call MPI_IRecv(f.dv(g.nf_y+1, g.first_x), 1, f_grey_zone, d, g.id, mp_cw, reqsd(3), this.ier)
	end if



	if (l > -1) then
		call MPI_ISend(f.d(g.first_y, g.ns_x), g.fstep*y, mp_dp, l, g.id, mp_cw, s_reqsl(1), this.ier);
		call MPI_ISend(f.du(g.first_y, g.ns_x), g.fstep*y, mp_dp, l, g.id, mp_cw, s_reqsl(2), this.ier);
		call MPI_ISend(f.dv(g.first_y, g.ns_x), g.fstep*y, mp_dp, l, g.id, mp_cw, s_reqsl(3), this.ier);
	end if
	if (r > -1) then
		call MPI_ISend(f.d(g.first_y, g.nf_x+1 - g.bstep), g.bstep*y, mp_dp, r, r, mp_cw, s_reqsr(1), this.ier);
		call MPI_ISend(f.du(g.first_y, g.nf_x+1 - g.bstep), g.bstep*y, mp_dp, r, r, mp_cw, s_reqsr(2), this.ier);
		call MPI_ISend(f.dv(g.first_y, g.nf_x+1 - g.bstep), g.bstep*y, mp_dp, r, r, mp_cw, s_reqsr(3), this.ier);
	end if

	if ( u > -1 ) then
		call MPI_ISend(f.d(g.ns_y, g.first_x), 1, f_grey_zone, u, u, mp_cw, s_reqsu(1), this.ier)
		call MPI_ISend(f.du(g.ns_y, g.first_x), 1, f_grey_zone, u, u, mp_cw, s_reqsu(2), this.ier)
		call MPI_ISend(f.dv(g.ns_y, g.first_x), 1, f_grey_zone, u, u, mp_cw, s_reqsu(3), this.ier)
	end if
	if ( d > -1 ) then
		call MPI_ISend(f.d(g.nf_y+1 - g.bstep, g.first_x), 1, b_grey_zone, d, d, mp_cw, s_reqsd(1), this.ier)
		call MPI_ISend(f.du(g.nf_y+1 - g.bstep, g.first_x), 1, b_grey_zone, d, d, mp_cw, s_reqsd(2), this.ier)
		call MPI_ISend(f.dv(g.nf_y+1 - g.bstep, g.first_x), 1, b_grey_zone, d, d, mp_cw, s_reqsd(3), this.ier)
	end if

! Wait RECV
	if ( l > -1 ) then
		call MPI_Waitall(3, reqsl, statsl, this.ier)
	end if
	if ( r > -1 ) then
		call MPI_Waitall(3, reqsr, statsr, this.ier)
	end if
	if ( u > -1 ) then
		call MPI_Waitall(3, reqsu, statsu, this.ier)
	end if
	if ( d > -1 ) then
		call MPI_Waitall(3, reqsd, statsd, this.ier)
	end if
! Wait SEND
	if ( l > -1 ) then
		call MPI_Waitall(3, s_reqsl, s_statsl, this.ier)
	end if
	if ( r > -1 ) then
		call MPI_Waitall(3, s_reqsr, s_statsr, this.ier)
	end if
	if ( u > -1 ) then
		call MPI_Waitall(3, s_reqsu, s_statsu, this.ier)
	end if
	if ( d > -1 ) then
		call MPI_Waitall(3, s_reqsd, s_statsd, this.ier)
	end if

end if

	End Subroutine


	Subroutine met_Print(this, f, g, t, name, Tmax, Wid, xid, yid, tid, ncid)

		Class(met) :: this
		Class(func) :: f
		Class(grid) :: g

		Integer x, y, i, j
		Integer status 
		Integer(4), Intent(In) :: t, Tmax, Wid, xid, yid, tid, ncid
		character(40), Intent(In) :: name

		Real(8) W_mass(g.ns_y:g.nf_y, g.ns_x:g.nf_x)

		do j = g.ns_x, g.nf_x
			do i = g.ns_y, g.nf_y
				W_mass(i, j) = f.d(i, j)
			end do
		end do


! 		if ( g.id == 0 ) then
		! ! указатель на первый элемент массива varval , число элементов
! 		  status = nf90_put_var (ncid, Wid, W_mass, (/ g.ns_y, g.ns_x, t/), (/ g.nf_y - g.ns_y + 1, g.nf_x - g.ns_x + 1, 1/) )
! 		end if

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