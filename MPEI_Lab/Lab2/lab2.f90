program lab2

	Use MPI

	Integer(4) lenght, first, last, id, np, ier, i, j, k, flag, p, q
! 	Lenght of massives, first and last index, process id, number of processes, error code
	Integer status(MPI_Status_size)
	Real(8) time_s, time_f
! 	Start time, finish time
	Real(8), Allocatable :: Y1(:), Y2(:), Y3(:), RESULT(:,:), A(:), B(:), C(:), D(:), E(:)
! 	Massives of exircise
	Real(8), Allocatable :: Temp1(:), Temp2(:), Temp3(:), Temp4(:), Temp5(:)


	lenght=16

	call MPI_Init(ier)
	call MPI_Comm_rank(MPI_Comm_World, id, ier)
	call MPI_Comm_size(MPI_Comm_World, np, ier)

	do j = 1, 4
		call MPI_Barrier(MPI_Comm_World, ier)
		time_s = MPI_Wtime()
		lenght=lenght*4


		first= id*lenght/np + 1 ! First index for process
		last= (id+1)*lenght/np

		Allocate(Y1(first:last), Y2(first:last), Y3(first:last), RESULT(1:lenght,1:3))
		Allocate(A(first:last), B(first:last), C(first:last), D(first:last),  E(first:last))
		Allocate(Temp1(first:last), Temp2(first:last), Temp3(first:last), Temp4(first:last), Temp5(first:last))

		do p = 1, 100, 1

			do i = first, last
				call random_number(A(i)); call random_number(B(i))
				call random_number(C(i)); call random_number(D(i))
				call random_number(E(i))
				A(i) = A(i)*10; B(i) = B(i)*10; C(i) = C(i)*10; D(i) = D(i)*10; E(i) = E(i)*10

				Temp1(i) = B(i) - A(i)
				Temp2(i) = C(i)*B(i)
				Temp3(i) = E(i) - A(i)/D(i)
				Temp4(i) = D(i)/A(i)
				Temp5(i) = A(i) + E(i) - C(i)
				Y2(i) = (Temp1(i)**2)*Temp2(i)
				Y1(i) = Temp1(i)*Temp2(i) + Temp4(i) - Y2(i)
				Y3(i) = Y2(i) + Temp3(i)

			end do

			if ( np > 1 ) then

				call MPI_Gather(Y1(first), last-first, MPI_DOUBLE_PRECISION, RESULT(1,1), last-first, MPI_DOUBLE_PRECISION, 0, MPI_Comm_World, ier)
				call MPI_Gather(Y2(first), last-first, MPI_DOUBLE_PRECISION, RESULT(1,2), last-first, MPI_DOUBLE_PRECISION, 0, MPI_Comm_World, ier)
				call MPI_Gather(Y3(first), last-first, MPI_DOUBLE_PRECISION, RESULT(1,3), last-first, MPI_DOUBLE_PRECISION, 0, MPI_Comm_World, ier)

			end if

		end do

		call MPI_Barrier(MPI_Comm_World, ier)

		if ( Allocated(Y1) ) then
			Deallocate(Y1); Deallocate(Y2); Deallocate(Y3); Deallocate(RESULT)
			Deallocate(A); Deallocate(B); Deallocate(C); Deallocate(D); Deallocate(E)
			Deallocate(Temp1); Deallocate(Temp2); Deallocate(Temp3);  Deallocate(Temp4);  Deallocate(Temp5)
		end if


		call MPI_Barrier(MPI_Comm_World, ier)
		time_f = MPI_Wtime()
		if ( id==0 ) then
			print *, time_f - time_s, lenght
		end if

	end do
	call MPI_Finalize(ier)

end

