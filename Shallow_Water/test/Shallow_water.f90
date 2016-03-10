program uravnenie

	Use modfunc, Only: func
	Use modnet, Only: grid
	Use method, Only: met
	Use Schema, Only: sch

!   Use netcdf

! 	Use SemiLagr, Only: Lagr
! 	Use Interp, Only: intp

	Implicit None

	include"mpif.h"

	integer x,y,t, ier, id,np, synchr, i, rc, eqvtype, iter_stpes, index, testid
	integer Xmax, Ymax, Tmax, bstep, fstep, timeset, casenumb, rcase, lcase, newfile
	integer Wid, xid, yid, tid

	character(40) name, sch_name
	character(1) str1
	character(2) str2
	real(8)  lengthX, lengthY, lengthT, cor, grav, height
	real(8), Allocatable :: time(:)
	real(8), Allocatable :: msgtime(:)
	real(8), Allocatable :: diftime(:)
	real(8), Allocatable :: msgdiftime(:)
	real(8), Allocatable :: kurant(:,:)
! 	real(8), Allocatable :: test(:,:,:)




	Type(grid) :: g
	Type(func) :: f, fprev, fnext
	Type(met) :: m
! 	Type(Lagr) :: l
	Type(sch) :: s
! 	Type(intp) :: inter

	integer status(MPI_STATUS_SIZE)  , ncid


	call MPI_Init(ier)
	call MPI_Comm_rank(MPI_COMM_WORLD,id,ier)
	call MPI_Comm_size(MPI_COMM_WORLD,np,ier)
! 	cw = MPI_COMM_WORLD



	lengthX=1;	lengthY=1;	lengthT=1.5
	Xmax = 100; Ymax = 100;	Tmax = 100;	synchr = 0
	bstep = 2;	fstep = 2
	rcase = 3;	timeset = 1;	eqvtype = 1;	casenumb = 1; iter_stpes = 2; lcase = 1
	newfile=0
! 	temp=1.35
	cor = 0
	grav = .2
	t = 0




	if ( id == 0 ) then

		open(9,file='/data4t/avolkov/Fortran/Shallow_Water/test/init.file')
			read(9, FMT="(5(I14), 4(e20.12))") Xmax, Ymax, Tmax, timeset, casenumb, grav, cor, height, lengthT
		close(9)

	end if


	call MPI_Bcast( Xmax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)
	call MPI_Bcast( Ymax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)
	call MPI_Bcast( Tmax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)
	call MPI_Bcast(timeset, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)
	call MPI_Bcast(casenumb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)
	call MPI_Bcast(grav, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
	call MPI_Bcast(cor, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
	call MPI_Bcast(height, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
	call MPI_Bcast(lengthT, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
	call mpi_barrier(MPI_COMM_WORLD, ier)


	Allocate(time(0:Tmax+1))
	Allocate(msgtime(1:2*Tmax))
	Allocate(diftime(1:Tmax))
	Allocate(msgdiftime(1:Tmax))



	if (casenumb == 1 ) then
		fstep = 1
		name = 'Linear.nc'; sch_name = 'Linear'
	elseif (casenumb == 5 ) then
		fstep = 4
		bstep = 4
		name = 'RungeK.nc'; sch_name = 'RungeKutta'
	end if

	call g.init (lengthX, lengthY, lengthT, Xmax, Ymax, Tmax, bstep, fstep, np, id)

	call f.init (g)
	call fprev.init (g)

	call m.init(g, status, ier, newfile)
	call s.init(casenumb, eqvtype, rcase, cor, grav, height, g)



	do x = f.ns_x, f.nf_x
		do y = f.ns_y, f.nf_y
			fprev.d(y, x) = height*exp(-((((10.0/Xmax)*(x-Xmax*0.5))**2)+(((10.0/Ymax)*(y-Ymax*0.5))**2))) + 0
			!f.d(y, x) = exp(-((((20.0/Xmax)*(x-Xmax/2))**2)+(((20.0/Ymax)*(y-Ymax/2))**2)))
! 			fprev.d(y,x) = id+1
! 		fprev.d(y, x) = id * 10
			fprev.du(y,x) = 0
			fprev.dv(y,x) = 0
		end do
	end do

! 	call m.to_print(fprev, g, 1, folder_name, 1, index)

! 	call m.BornParam(fprev, g)


! 		status = nf90_create (path = trim('/data4t/avolkov/Fortran/Shallow_Water/datFiles/'//name), &
! cmode = IOR(NF90_NETCDF4,IOR(NF90_MPIIO,NF90_CLOBBER)), comm = MPI_COMM_WORLD, info = MPI_INFO_NULL, ncid = ncid)

! 		status = nf90_def_dim (ncid, "x", g.StepsX, xid)
! 		status = nf90_def_dim (ncid, "y", g.StepsY, yid)
! 		status = nf90_def_dim (ncid, "t", Tmax, tid)

! 		status = nf90_def_var (ncid, "water", NF90_REAL, (/ yid, xid, tid/), Wid)
! 		status = nf90_enddef  (ncid)
! 		call m.to_print(fprev, g, t, name, Tmax, Wid, xid, yid, tid, ncid)

	index = 1

	call mpi_barrier(MPI_COMM_WORLD, ier)
	if ( id == 0 ) then
		time(0) = MPI_Wtime()
	end if

	do t=1,Tmax

		SELECT CASE (casenumb)
		CASE(0)

		CASE (1)
			call s.linear(f.d, fprev.d, f.du, fprev.du, f.dv, fprev.dv, g)
			call m.Message(f, g)!;call m.Message(f.du, g);call m.Message(f.dv, g)

		CASE(5)
			call s.RungeKutta(f, fprev, g, m)
			msgtime(t) = MPI_Wtime()
			call m.Message(f, g)!;call m.Message(f.du, g);call m.Message(f.dv, g)
			msgtime(2*t) = MPI_Wtime()
		END SELECT 


		do x = g.ns_x, g.nf_x
			if ( g.ns_y == 1 ) then
				f.dv(g.ns_y,x) = -f.dv(g.ns_y,x)
				f.du(g.ns_y,x) = -f.du(g.ns_y,x)
			end if

			if ( g.nf_y == g.StepsY ) then
				f.dv(g.nf_y,x) = -f.dv(g.nf_y,x)
				f.du(g.nf_y,x) = -f.du(g.nf_y,x)
			end if
		end do

		do y = g.ns_y, g.nf_y
			if ( g.ns_x == 1 ) then
				f.dv(y,g.ns_x) = -f.dv(y,g.ns_x)
				f.du(y,g.ns_x) = -f.du(y,g.ns_x)
			end if

			if ( g.nf_x == g.StepsX ) then
				f.dv(y,g.nf_x) = -f.dv(y,g.nf_x)
				f.du(y,g.nf_x) = -f.du(y,g.nf_x)
			end if
		end do

		call fprev.eq(f)
! 		call m.to_print(fprev, g, t, name, Tmax, Wid, xid, yid, tid, ncid)
! 		end if

		call mpi_barrier(MPI_COMM_WORLD, ier)
		if ( id == 0 ) then
			time(t) = MPI_Wtime()
			diftime(t) = time(t) - time(t-1)
			msgdiftime(t) = msgtime(2*t) - msgtime(t)
		end if

	end do




! 	print *,11,id

	call mpi_barrier(MPI_COMM_WORLD, ier)

! 	status = nf90_close (ncid)

	if ( g.id == 0 ) then
		print *,"done"
! 		print *,"dx", g.dx, "dy", g.dy, "dt", g.dt
		print *, "Your file is: "//name
	end if

	call f.deinit()
	call fprev.deinit()
	call s.deinit()
! 	call l.deinit()

	call mpi_barrier(MPI_COMM_WORLD, ier)
	if ( id == 0 ) then
! 		do t = 1,Tmax
! 			print *, "cycle time = ", diftime(t)
! 		end do
		time(Tmax+1) = MPI_Wtime()
		print *, "time = ", time(Tmax+1) - time(0)
		print *, "max cycle time = ", MAXVAL(diftime)
		print *, "min cycle time = ", MINVAL(diftime)
		print *, "med cycle time = ", (time(Tmax+1) - time(0))/Tmax
		print *, "Msg all time = ", SUM(msgdiftime)
	end if



	call MPI_FINALIZE(rc)

end
