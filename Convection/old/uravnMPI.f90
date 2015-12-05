program uravnenie

	Use MPI
	Use modfunc, Only: func
	Use modnet, Only: grid
	Use method, Only: met

	Implicit None


	integer x,t, ier, id,np, synchr, i, rc, eqvtype
	integer Xmax, Tmax, param, bstep, fstep, timeset, casenumb
	character(20) name
	real(8) u, lengthX, lengthT

	Type(grid) :: g
	Type(func) :: f, fprev, fnext
	Type(met) :: m

	integer status(MPI_STATUS_SIZE)  


	call MPI_Init(ier)

	call MPI_Comm_rank(MPI_COMM_WORLD,id,ier)
	call MPI_Comm_size(MPI_COMM_WORLD,np,ier)
	

	lengthX=10
	lengthT=10
	Xmax = 100
	Tmax = 100
	u=.2
	bstep = 1
	fstep = 0
	synchr = 0
	casenumb = 1
	timeset = 1
	eqvtype = 0

	if ( casenumb == 1 ) then
		fstep = 0
		name = 'Angle.dat'
	elseif ( casenumb == 2 ) then
		fstep = 1
		name = 'Cross.dat'
	elseif ( casenumb == 3 ) then
		fstep = 1
		name = 'Lax.dat'
	elseif ( casenumb == 4 ) then
		fstep = 1
		name = 'LaxVen.dat'
	end if

	if ( eqvtype == 0 ) then
		print *, "Linear equation"
	elseif ( eqvtype == 1 ) then
		print *, "Burgers"
	end if

		call g.init (lengthX, lengthT, Xmax, Tmax, bstep, fstep, u, np, id)
		call m.init(casenumb)
		call f.init (g.ns-g.bstep,g.nf+g.fstep)
		call fprev.init (g.ns-g.bstep,g.nf+g.fstep)

		if ( casenumb == 2 .OR. casenumb == 4 ) then
			call fnext.init(g.ns-bstep,g.nf+fstep)
		end if



		do x = g.ns-bstep, g.nf+fstep
			fprev.d(x) = exp(-0.02*(x-Xmax/2)**2)
		end do


		write(*,'(2x,a,i4,3x,a,i4)') "myid =",id,"np=",np
		print *, g.ns, g.nf 




	SELECT CASE (casenumb)
		CASE (1)

			if (id == 0 ) then
				print *, "Angle scheme!"
			end if

			do t=1,Tmax-1

				call m.Message(fprev, g, status, ier)
				call m.angle(f, fprev, g, eqvtype)
				call m.to_print(f, g, t, status, ier, name, timeset)

			end do


		CASE (2)

			if (id == 0 ) then
				print *, "Cross scheme!"
			end if

			call m.angle(f, fprev, g, eqvtype)


			do t=2,Tmax-1

				call m.Message(f, g, status, ier)
				call m.cross(f, fprev, fnext, g, eqvtype)
				call m.to_print(f, g, t, status, ier, name, timeset)

			end do


		CASE (3)

			if ( id == 0 ) then
				print *, "Lax scheme!"
			end if

			do t=1,Tmax-1


				call m.Lax(f, fprev, g)
				!call m.Message(f, g, status, ier)
				!call fprev.eq(f)
				call m.to_print(f, g, t, status, ier, name, timeset)

			end do


		CASE (4)

			if ( id == 0 ) then
				print *, "LaxVendroff scheme!"
			end if

			do t=1,Tmax-1

				if ( mod(t,2)==0  ) then

					call m.Message(f, g, status, ier)

				else

					call m.Message(fprev, g, status, ier)

				end if
				
				call m.LaxVen(f, fprev, fnext, g, t)
				call m.to_print(fprev, g, t, status, ier, name, timeset)

			end do

	END SELECT 



	print *,11,id

	call mpi_barrier(MPI_COMM_WORLD, ier)
	print *,"done"
	

	call MPI_FINALIZE(rc)
	 
end

