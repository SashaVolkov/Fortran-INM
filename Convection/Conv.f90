program uravnenie

	Use MPI
	Use modfunc, Only: func
	Use modnet, Only: grid
	Use method, Only: met
	Use Schema, Only: sch
	Use SemiLagr, Only: Lagr
	Use Interp, Only: intp

	Implicit None


	integer x,t, ier, id,np, synchr, i, rc, eqvtype, iter_stpes, Res_first, Res_last
	integer Xmax, Tmax, param, bstep, fstep, timeset, casenumb, rcase, lcase, newfile
	character(40) name, sch_name
	real(8) u, lengthX, lengthT, pi


	Type(grid) :: g
	Type(func) :: f, fprev, fnext
	Type(met) :: m
	Type(Lagr) :: l
	Type(sch) :: s
	Type(intp) :: inter

	integer status(MPI_STATUS_SIZE)  


	call MPI_Init(ier)
	call MPI_Comm_rank(MPI_COMM_WORLD,id,ier)
	call MPI_Comm_size(MPI_COMM_WORLD,np,ier)

	lengthX=1;	lengthT=1
	Xmax = 10;	Tmax = 10;	synchr = 0
	bstep = 1;	fstep = 0
	u=2.5
	rcase = 3;	timeset = 0;	eqvtype = 0;	casenumb = 0; iter_stpes = 2
	newfile=0
! 	temp=1.35



	if ( id == 0 ) then

		open(9,file='/home/sasha/Fortran/Convection/init.file')
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
		u = .2
	end if

	if ( casenumb == 0 ) then
		name = 'test.dat';
		bstep = 5
		fstep = 5
	end if

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
		fstep = 3+int(abs(u))
		bstep = 3+int(abs(u))
		name = 'Lagr.dat'; sch_name = 'Semi-Lagrangian'
	end if

	call g.init (lengthX, lengthT, Xmax, Tmax, bstep, fstep, u, np, id)
	call f.init (g.ns-g.bstep,g.nf+g.fstep)
	call fprev.init (g.ns-g.bstep,g.nf+g.fstep)
	call m.init(g, status, ier, newfile)
	call s.init(casenumb, eqvtype, rcase, g)
	call l.init(lcase, iter_stpes, eqvtype, g)



	if ( eqvtype == 1 .AND. g.id == 0) then
		print *, "Linear equation"
	elseif ( eqvtype == 2 .AND. g.id == 0) then
		print *, "Burgers"
	end if


	if ( (u*g.dt/g.dx) > 1 .and. casenumb /=6 .and. casenumb /=0) then
		print *, 'CFL condition is too bad! Change quantity of steps, or change the velocity!'
		call MPI_FINALIZE(rc)
		STOP
	end if


	if ( casenumb == 2 .OR. casenumb == 4) then
		call fnext.init(g.ns-bstep,g.nf+fstep)
	end if


	do x = g.ns-bstep, g.nf+fstep
		fprev.d(x) = exp(-((20.0/Xmax)*(x-Xmax/2))**2)

! 		fprev.d(x) = id + 1
	end do

	call m.to_print(fprev, g, 1, name, 1)



	call m.BornParam(fprev, g)


	SELECT CASE (casenumb)
	CASE(0)

	CASE (1)

		if (id == 0 ) then
			print *, "Angle scheme!"
		end if

		do t=1,Tmax-1

			call m.Message(fprev.d, g)
			call s.angle(f.d, fprev.d, g)
			call fprev.eq(f)
			if(t == timeset)then
				call m.to_print(f, g, t, name, timeset)
			end if

		end do


	CASE (2)

		if (id == 0 ) then
			print *, "Cross scheme!"
		end if

		call s.angle(f.d, fprev.d, g)
		call fprev.eq(f)


		do t=2,Tmax-1
			if ( g.np > 1 ) then
				call m.Message(f.d, g)
			end if
			call s.cross(f.d, fprev.d, fnext.d, g)
			call fprev.eq(f)
			call f.eq(fnext)
			if(t == timeset)then
				call m.to_print(f, g, t, name, timeset)
			end if
		end do


	CASE (3)

		if ( id == 0 ) then
			print *, "Lax scheme!"
		end if

		do t=1,Tmax-1
			call s.Lax(f.d, fprev.d, g)
			call m.Message(f.d, g)
			call fprev.eq(f)
			if(t == timeset)then
				call m.to_print(f, g, t, name, timeset)
			end if
		end do


	CASE (4)

		if ( id == 0 ) then
			print *, "LaxVendroff scheme!"
		end if

		do t=1,Tmax-1
			call s.LaxVen(f, fprev, fnext, g, m, t)
			if(t == timeset)then
				call m.to_print(fnext, g, t, name, timeset)
			end if
		end do


	CASE (5)

		if ( id == 0 ) then
			print *, "Runge Kutta"
		end if

		do t=1,Tmax-1
			call m.Message(fprev.d, g)
			call s.RungeKutta(f, fprev, g, m)
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
				call l.SemiLagr(fprev.d, f.d, g, inter)
				call fprev.eq(f)
				if(t == timeset)then
					call m.to_print(fprev, g, t, name, timeset)
				end if
			end do




	END SELECT 

	if ( casenumb /= 0 ) call m.SchemeParam(f, g, sch_name)

	print *,11,id

	call mpi_barrier(MPI_COMM_WORLD, ier)

	if ( g.id == 0 ) then
		print *,"done"
		print *, "Your file is: "//name
	end if

	call f.deinit()
	call fprev.deinit()
	call s.deinit()
	call l.deinit()

	call MPI_FINALIZE(rc)

end
