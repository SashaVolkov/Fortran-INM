program uravnenie

	Use MPI
	Use modfunc, Only: func
	Use modnet, Only: grid
	Use method, Only: met
	Use Schema, Only: sch
! 	Use SemiLagr, Only: Lagr
! 	Use Interp, Only: intp

	Implicit None


	integer x,y,t, ier, id,np, synchr, i, rc, eqvtype, iter_stpes, index
	integer bstep, fstep, timeset, casenumb, rcase, lcase, newfile
	character(6) folder_name, sch_name
	real(8)  lengthX, lengthY, lengthT
	real(8), Allocatable :: kurant(:,:)



	Type initial
		integer :: Xmax, Ymax, Tmax, timeset, casenumb
		real(8) :: cor, grav, height
	End Type initial


	Type(grid) :: g
	Type(func) :: f, fprev, fnext
	Type(met) :: m
	Type(sch) :: s
	Type(initial) :: init

	integer status(MPI_STATUS_SIZE)


	call MPI_Init(ier)
	call MPI_Comm_rank(MPI_COMM_WORLD,id,ier)
	call MPI_Comm_size(MPI_COMM_WORLD,np,ier)

	lengthX=1; lengthY=1; lengthT=2; index=0
	synchr = 0; bstep = 1; fstep = 0
	rcase = 3; timeset = 1; eqvtype = 1; casenumb = 1;
	iter_stpes = 2; lcase = 1; newfile=0


!initialization
! 	if ( id == 0 ) then
		open(9,file='/home/sasha/Fortran/Shallow_Water/test1/init.file')
! 	end if
			read(9, FMT="(5(I14), 3(e20.12))") init.Xmax, init.Ymax, init.Tmax,&
			init.timeset, init.casenumb, init.grav, init.cor, init.height
	if ( id == 0 ) then
		close(9)
	end if
! 	call MPI_Bcast( init, 12, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)
! 	call mpi_barrier(MPI_COMM_WORLD, ier)
!end of initialization


! if ( id == 1 ) then
! 	print *, init.Xmax, init.Ymax, init.Tmax,&
! 			init.timeset, init.casenumb, init.grav, init.cor, init.height
! end if

	if ( init.casenumb == 1 ) then
		fstep = 1
		folder_name = 'Linear'; sch_name = 'Linear'
	elseif ( init.casenumb == 5 ) then
		fstep = 1
		folder_name = 'RungeK'; sch_name = 'RungeKutta'
	end if




	call g.init (lengthX, lengthY, lengthT, init.Xmax, init.Ymax, init.Tmax, bstep, fstep, np, id)
	call f.init (g)
	call fprev.init (g)
	call m.init(g, status, ier, newfile)
	call s.init(casenumb, eqvtype, rcase, init.cor, init.grav, init.height, g)


! first moment state
	do x = f.ns_x, f.nf_x
		do y = f.ns_x, f.nf_x
			fprev.d(y, x) = exp(-((((10.0/init.Xmax)*(x-init.Xmax*0.5))**2)+(((10.0/init.Ymax)*(y-init.Ymax*0.5))**2))) + 0
			!f.d(y, x) = exp(-((((20.0/Xmax)*(x-Xmax/2))**2)+(((20.0/Ymax)*(y-Ymax/2))**2)))
			fprev.du(y,x) = 0
			fprev.dv(y,x) = 0
		end do
	end do

	call m.to_print(fprev, g, 1, folder_name, 1, index)

! 	call m.BornParam(fprev, g)
	index = 1

	do t=1,init.Tmax-1

		call m.Message(fprev.d, g)
		call m.Message(fprev.du, g)
		call m.Message(fprev.dv, g)

			SELECT CASE (casenumb)
			CASE(0)

			CASE (1)
				call s.linear(f.d, fprev.d, f.du, fprev.du, f.dv, fprev.dv, g)
			CASE(5)
				call s.RungeKutta(f, fprev, g, m)
			END SELECT 


			do x = f.ns_x, f.nf_x
				f.dv(g.ns_y,x) = -f.dv(g.ns_y,x)
				f.dv(f.ns_y,x) = -f.dv(f.ns_y,x)
				f.dv(g.nf_y,x) = -f.dv(g.nf_y,x)
				f.dv(f.nf_y,x) = -f.dv(f.nf_y,x)
			end do

			do y = f.ns_y, f.nf_y
				f.du(y,g.ns_x) = -f.du(y,g.ns_x)
				f.du(y,f.ns_x) = -f.du(y,f.ns_x)
				f.du(y,g.nf_x) = -f.du(y,g.nf_x)
				f.du(y,f.nf_x) = -f.du(y,f.nf_x)
			end do



		call fprev.eq(f)


		if(t == index*init.timeset)then
			if ( id == 0 ) then
				print *, t
			end if

			call m.to_print(f, g, t, folder_name, init.timeset, index)
			index = index + 1
			if ( index == 10 ) then
				index = 0
			end if
		end if

	end do


! 	print *,11,id

	call mpi_barrier(MPI_COMM_WORLD, ier)

	if ( g.id == 0 ) then
		print *,"done"
! 		print *,"dx", g.dx, "dy", g.dy, "dt", g.dt
		print *, "Your folder is: "//folder_name
	end if

	call f.deinit()
	call fprev.deinit()
	call s.deinit()
! 	call l.deinit()

	call MPI_FINALIZE(rc)

end
