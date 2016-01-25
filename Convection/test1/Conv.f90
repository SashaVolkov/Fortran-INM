program uravnenie

	Use MPI
	Use modfunc, Only: func
	Use modnet, Only: grid
	Use method, Only: met
	Use Schema, Only: sch
! 	Use SemiLagr, Only: Lagr
! 	Use Interp, Only: intp

	Implicit None


	integer x,y,t, ier, id,np, synchr, i, rc, eqvtype, iter_stpes
	integer Xmax, Ymax, Tmax, bstep, fstep, timeset, casenumb, rcase, lcase, newfile

	character(40) name, sch_name
	real(8) u_x, u_y, lengthX, lengthY, lengthT, pi


	Type(grid) :: g
	Type(func) :: f, fprev, fnext
	Type(met) :: m
	Type(Lagr) :: l
	Type(sch) :: s
! 	Type(intp) :: inter

	integer status(MPI_STATUS_SIZE)  


	call MPI_Init(ier)
	call MPI_Comm_rank(MPI_COMM_WORLD,id,ier)
	call MPI_Comm_size(MPI_COMM_WORLD,np,ier)

	lengthX=1;	lengthY=1;	lengthT=1
	Xmax = 100; Ymax = 100;	Tmax = 100;	synchr = 0
	bstep = 1;	fstep = 0
	u_x=.7; u_y = .7
	rcase = 3;	timeset = 1;	eqvtype = 1;	casenumb = 1; iter_stpes = 2; lcase = 1
	newfile=0
! 	temp=1.35



	if ( id == 0 ) then

		open(9,file='/home/sasha/Fortran/Convection/test1/init.file')
			read(9, FMT="(6(I14))") Xmax, Tmax, timeset, eqvtype, casenumb, lcase
		close(9)

	end if


	call MPI_Bcast(casenumb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)
	call MPI_Bcast( Xmax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)
	call MPI_Bcast( Tmax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)
	call MPI_Bcast(lcase, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)
	call MPI_Bcast(timeset, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)
	call MPI_Bcast(eqvtype, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)
	call mpi_barrier(MPI_COMM_WORLD, ier)



	if(casenumb > 0 .and. casenumb < 6) then
		u_x = .7
		u_y = .7
	end if

! 	if ( casenumb == 0 ) then
! 		name = 'test.dat';
! 		bstep = 5
! 		fstep = 5
! 	end if

	if ( casenumb == 1 ) then
		fstep = 0
		name = 'Angle.dat'; sch_name = 'Angle'
	elseif ( casenumb == 2 ) then
		fstep = 1
		name = 'Cross.dat'; sch_name = 'Cross'
	elseif ( casenumb == 3 ) then
		fstep = 1
		name = 'Lax.dat'; sch_name = 'Lax'
	elseif ( casenumb == 4 ) then
		fstep = 1
		name = 'LaxVen.dat'; sch_name = 'LaxVendroff'
	elseif ( casenumb == 5 ) then
	if ( rcase == 1 ) then
		fstep = 0
	elseif ( rcase == 2 ) then
		fstep = 1
	elseif ( rcase == 3 ) then
		fstep = 2
		bstep = 2
	end if
		name = 'RK.dat'; sch_name = 'RungeKutta'
	elseif ( casenumb == 6 ) then
		fstep = 3+int(abs(u_x))
		bstep = 3+int(abs(u_x))
		name = 'Lagr.dat'; sch_name = 'Semi-Lagrangian'
	end if

	call g.init (lengthX, lengthY, lengthT, Xmax, Ymax, Tmax, bstep, fstep, u_x, u_y, np, id)

	call f.init (g)
	call fprev.init (g)

	call m.init(g, status, ier, newfile)
	call s.init(casenumb, eqvtype, rcase, g)
	call l.init(lcase, iter_stpes, eqvtype, g)



	if ( eqvtype == 1 .AND. g.id == 0) then
		print *, "Linear equation"
	elseif ( eqvtype == 2 .AND. g.id == 0) then
		print *, "Burgers"
	end if


	if ( (u_x*g.dt/g.dx) > 1 .and. casenumb /=6 .and. casenumb /=0) then
		print *, 'CFL condition is too bad! Change quantity of steps, or change the velocity!'
		call MPI_FINALIZE(rc)
		STOP
	end if


! 	if ( casenumb == 2 .OR. casenumb == 4) then
! 		call fnext.init(g)
! 	end if


	do x = g.ns_x-bstep, g.nf_x+fstep
		do y = 1, g.StepsY
			fprev.d(y, x) = exp(-((((20.0/Xmax)*(x-Xmax/2))**2)+(((20.0/Ymax)*(y-Ymax/2))**2)))

! 					fprev.d(y,x) = id + 1

		end do
	end do

	call m.to_print(fprev, g, 1, name, 1)



! 	call m.BornParam(fprev, g)


	SELECT CASE (casenumb)
	CASE(0)

	CASE (1)

		if (id == 0 ) then
			print *, "Angle scheme!"
		end if

		do t=1,Tmax-1

			call m.Message(fprev.d, g)
			do x = g.ns_x-bstep, g.nf_x+fstep
				fprev.d(1, x) = fprev.d(g.StepsY,x)
			end do
			call s.angle(f.d, fprev.d, g)
			call fprev.eq(f)
			if(t == timeset)then
				call m.to_print(f, g, t, name, timeset)
			end if

		end do



		CASE (6)

			if ( id == 0 ) then
				print *, "Semi-Lagrangian"
			end if


			do t=1, Tmax-1
				call m.Message(fprev.d, g)

				do x = g.ns_x-bstep, g.nf_x+fstep
					fprev.d(1, x) = fprev.d(g.StepsY,x)
				end do


				
				call l.SemiLagr(fprev.d, f.d, g, inter)
				call fprev.eq(f)
				if(t == timeset)then
					call m.to_print(fprev, g, t, name, timeset)
				end if
			end do



	END SELECT 

! 	if ( casenumb /= 0 ) call m.SchemeParam(f, g, sch_name)

	print *,11,id

	call mpi_barrier(MPI_COMM_WORLD, ier)

	if ( g.id == 0 ) then
		print *,"done"
		print *, "Your file is: "//name
	end if

	call f.deinit()
	call fprev.deinit()
	call s.deinit()
! 	call l.deinit()

	call MPI_FINALIZE(rc)

end
