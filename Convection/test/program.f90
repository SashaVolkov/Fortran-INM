program test
	Use MPI
	Use modnet, Only: grid
	Use method, Only: met
	Use Interp, Only: intp



		Implicit None


	integer x,t, ier, id,np, synchr, i, rc, iter_stpes, Res_first, Res_last
	integer Xmax, Tmax, param, bstep, fstep, timeset, casenumb, rcase, lcase, newfile
	character(40) name, sch_name
	real(8) u, lengthX, lengthT, pi
! 	real(4) :: temp
	Real(8), Allocatable :: RESULT(:,:)
	Real(8), Allocatable :: TEMP(:,:)
	Real(8), Allocatable :: TEMP1(:,:)

	Type(grid) :: g
	Type(met) :: m
	Type(intp) :: inter

	integer status(MPI_STATUS_SIZE)  


	call MPI_Init(ier)
	call MPI_Comm_rank(MPI_COMM_WORLD,id,ier)
	call MPI_Comm_size(MPI_COMM_WORLD,np,ier)

	lengthX=1;	lengthT=1
	Xmax = 10;	Tmax = 10;	synchr = 0
	bstep = 1;	fstep = 0
	u=2.5
	rcase = 3;	timeset = 0;	casenumb = 0; iter_stpes = 2
	newfile=0
! 	temp=1.35


	call mpi_barrier(MPI_COMM_WORLD, ier)


	if ( casenumb == 0 ) then
		name = 'test.dat';
		bstep = 5
		fstep = 5
	end if


	call g.init (lengthX, lengthT, Xmax, Tmax, bstep, fstep, u, np, id)
	call m.init(g, status, ier, newfile)

	Res_first = 1
	Res_last = Xmax*4

	Allocate(RESULT(Res_first-12:Res_last+12, Res_first-12:Res_last+12))
	Allocate(TEMP(g.ns-g.bstep:g.nf+g.fstep, g.ns-g.bstep:g.nf+g.fstep))
	Allocate(TEMP1(g.ns-g.bstep:g.nf+g.fstep, Res_first-12:Res_last+12))


	do x = g.ns-bstep, g.nf+fstep
		do i = g.ns-bstep, g.nf+fstep
			TEMP(x,i) = x + i
		end do
	end do


	do x=g.ns-g.bstep, g.nf+g.fstep
		do i = Res_first, Res_last
			call inter.Lintp(real(i/4.0,8), TEMP1(x,i), TEMP(x,:), g.ns-g.bstep, g.nf+g.fstep)
		end do
	end do


	do i = Res_first, Res_last
			do x = Res_first, Res_last
			call inter.Lintp(real(x/4.0,8), RESULT(x,i), TEMP1(:,i), g.ns-g.bstep, g.nf+g.fstep)
		end do
	end do



	open(14,file='/home/sasha/Fortran/Convection/datFiles/interp.dat')
	do i = 1, Res_last
		write(14,'(3(3x,e20.12))') (Real(x-1,8)*g.dx/4.0, Real(i-1,8)*g.dx/4.0, RESULT(x,i), x=1,Res_last)
	end do
	close(14)

	open(14,file='/home/sasha/Fortran/Convection/datFiles/interp1.dat')
	do i = 1, Xmax
		write(14,'(3(3x,e20.12))') (Real(x-1,8)*g.dx, Real(i-1,8)*g.dx, TEMP(x,i), x=1,Xmax)
	end do
	close(14)

	call MPI_Finalize(ier)

end