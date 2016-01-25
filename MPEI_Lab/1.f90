program FirstLab

	Use MPI
	IMPLICIT NONE


	Real(8) S1, S2, time, time_start, time_finish
	Real(8), Allocatable :: A(:)
	Real(8), Allocatable :: B(:)
	Real(8), Allocatable :: C(:)
	Real(8), Allocatable :: Y(:)
	integer Size_All, first, last, j, i, param
	integer  id, np, ier, rc
	integer status(MPI_STATUS_SIZE)  

	call MPI_Init(ier)
	call MPI_Comm_rank(MPI_COMM_WORLD,id,ier)
	call MPI_Comm_size(MPI_COMM_WORLD,np,ier)

	time_start = MPI_WTIME()

	Size_All = 8
	first = id*Size_All + 1
	last = (id + 1)*Size_All

	S1 = 33.6
	S2 = 45.2
	param = 1000000

	Allocate(A(first:last))
	Allocate(B(first:last))
	Allocate(C(first:last))
	Allocate(Y(first:last))

	do i = first, last
		
		A(i) = exp((i-Size_All)*1.0)
		B(i) = sin((i)*1.0)
		C(i) = asin((i**2)*1.0)

	end do

	do j = 1, param
		do i = first, last

			Y(i) = (A(i) + S1*B(i))*C(i) + S2
			
		end do
	end do


	if ( id == 0 ) then 
		param = 1
	    open(14,file='result')
	    write(14,'(2(3x,e20.12))') (Real(8)*i, Y(i), i=first,last)
	    close(14)

	    if ( np > 1 ) then
	     	call MPI_Send(param, 1, MPI_INTEGER, id+1, id+1, MPI_COMM_WORLD, ier)
	     end if 

    else
	    call MPI_Recv(param, 1, MPI_INTEGER, id-1, id, MPI_COMM_WORLD, status, ier)

	    open(14,file='result', position="append")
	    write(14,'(2(3x,e20.12))') (Real(8)*i, Y(i), i=first,last)
	    close(14)

	end if

	time_finish = MPI_WTIME()
	time = time_finish - time_start

	call MPI_Barrier(MPI_COMM_WORLD, ier)

	if ( id == 0 ) then
		
		print *, "Size = ", Size_All
		print *, "Time = ", time
		
	end if

	call MPI_FINALIZE(rc)


end