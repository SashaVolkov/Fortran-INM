program uravnenie

	Use MPI
	Use modfunc, Only: func
	Use modnet, Only: grid
	Use method, Only: met
	Use Schema, Only: sch

  Use netcdf

! 	Use SemiLagr, Only: Lagr
! 	Use Interp, Only: intp

	Implicit None

	integer x,y,t, ier, id,np, synchr, i, rc, eqvtype, iter_stpes, index, testid
	integer Xmax, Ymax, Tmax, bstep, fstep, timeset, casenumb, rcase, lcase, newfile

	character(40) name, sch_name
	character(1) str1
	character(2) str2
	real(8)  lengthX, lengthY, lengthT, cor, grav, height
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

	lengthX=1;	lengthY=1;	lengthT=2
	Xmax = 100; Ymax = 100;	Tmax = 100;	synchr = 0
	bstep = 1;	fstep = 0
	rcase = 3;	timeset = 1;	eqvtype = 1;	casenumb = 1; iter_stpes = 2; lcase = 1
	newfile=0
! 	temp=1.35
	cor = 0
	grav = .2




	if ( id == 0 ) then

		open(9,file='/home/sasha/Fortran/Shallow_Water/test/init.file')
			read(9, FMT="(5(I14), 3(e20.12))") Xmax, Ymax, Tmax, timeset, casenumb, grav, cor, height
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
	call mpi_barrier(MPI_COMM_WORLD, ier)


! 	Allocate(test(1:Xmax, 1:Ymax, 1:Tmax))

	if (casenumb == 1 ) then
		fstep = 1
		name = 'Linear.dat'; sch_name = 'Linear'
	elseif (casenumb == 5 ) then
		fstep = 1
		name = 'RungeK.dat'; sch_name = 'RungeKutta'
	end if

	call g.init (lengthX, lengthY, lengthT, Xmax, Ymax, Tmax, bstep, fstep, np, id)

	call f.init (g)
	call fprev.init (g)

	call m.init(g, status, ier, newfile)
	call s.init(casenumb, eqvtype, rcase, cor, grav, height, g)



	do x = f.ns_x, f.nf_x
		do y = f.ns_y, f.nf_y
			fprev.d(y, x) = height*exp(-((((20.0/Xmax)*(x-Xmax*0.5))**2)+(((20.0/Ymax)*(y-Ymax*0.5))**2))) + 0
			!f.d(y, x) = exp(-((((20.0/Xmax)*(x-Xmax/2))**2)+(((20.0/Ymax)*(y-Ymax/2))**2)))
			fprev.du(y,x) = 0
			fprev.dv(y,x) = 0
		end do
	end do

! 	call m.to_print(fprev, g, 1, folder_name, 1, index)

! 	call m.BornParam(fprev, g)


		status = nf90_create (path = trim("filename.nc"), &
cmode = IOR(NF90_NETCDF4,IOR(NF90_MPIIO,NF90_CLOBBER)), comm = MPI_COMM_WORLD, info = MPI_INFO_NULL, ncid = ncid)

	index = 1

	do t=1,Tmax

! 		call m.Message(fprev.d, g)
! 		call m.Message(fprev.du, g)
! 		call m.Message(fprev.dv, g)

			SELECT CASE (casenumb)
			CASE(0)

			CASE (1)
				call s.linear(f.d, fprev.d, f.du, fprev.du, f.dv, fprev.dv, g)
			CASE(5)
				call s.RungeKutta(f, fprev, g, m)
			END SELECT 


			do x = g.ns_x, g.nf_x
				f.dv(g.ns_y,x) = f.dv(f.nf_y,x)
				f.du(g.ns_y,x) = f.du(f.nf_y,x)
				f.d(g.ns_y,x) = f.d(f.nf_y,x)

				f.dv(f.ns_y,x) = f.dv(g.nf_y,x)
				f.du(f.ns_y,x) = f.du(g.nf_y,x)
				f.d(f.ns_y,x) = f.d(g.nf_y,x)
			end do


		call fprev.eq(f)
! 		if(t == index*timeset)then
				call m.to_print(fprev, g, t, name, Tmax)
! 		end if

	end do

! 	if ( casenumb /= 0 ) call m.SchemeParam(f, g, sch_name)

	print *,11,id

	call mpi_barrier(MPI_COMM_WORLD, ier)

	status = nf90_close (ncid)

	if ( g.id == 0 ) then
		print *,"done"
! 		print *,"dx", g.dx, "dy", g.dy, "dt", g.dt
		print *, "Your file is: "//name
	end if

	call f.deinit()
	call fprev.deinit()
	call s.deinit()
! 	call l.deinit()

	call MPI_FINALIZE(rc)

end
