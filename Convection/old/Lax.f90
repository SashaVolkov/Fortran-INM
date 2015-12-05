program Lax

  Use MPI

  Implicit None


  real(8) lengthX, lengthT, u, dx, dt, G
  integer x,t, ier, id,np, synchr, i, rc
  integer Xmax, Tmax, param, Xsize, first, last
  real(8), Allocatable :: func(:,:)

  integer status(MPI_STATUS_SIZE)  

  call MPI_Init(ier)

  call MPI_Comm_rank(MPI_COMM_WORLD,id,ier)
  call MPI_Comm_size(MPI_COMM_WORLD,np,ier)
  

  lengthX=10
  lengthT=10
  Xmax = 100
  Tmax = 100
  dx = lengthX/Xmax; dt = lengthT/Tmax
  u=.2
 
  synchr = 0
  Xsize = Xmax/np
  first = id*Xsize+1
  last = (id+1)*Xsize

  G=u*(dt/dx)

  if ( dt<dx/u*dx+0.005 ) then
    if ( dt>dx/u*dx-0.005 ) then
      write(*,*) "BULLSHIT!"      
    end if
  end if


  write(*,'(2x,a,i4,3x,a,i4,3x,a,f14.7)') "myid =",id,"np=",np,"dx=",dx


    Allocate(func(first-1:last+1, Tmax))

    do x = first-1, last+1
      func(x,1) = exp(-0.02*(x-Xmax/2)**2)
    end do
  

!    do x= first,last+1
!      func(x,2)=func(x,1)-u*(dt/dx)*(func(x,1)-func((x-1),1))
!    end do  


  print *,id



  do t=1,Tmax-1


    if ( id==np-1 ) then
      call MPI_Send(func( last, t), 1, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ier)
    else
      call MPI_Send(func( last, t), 1, MPI_DOUBLE_PRECISION, id+1, id+1, MPI_COMM_WORLD, ier)
    end if

    if ( id==0 ) then
      call MPI_Recv(func( first-1, t), 1, MPI_DOUBLE_PRECISION, np-1, id, MPI_COMM_WORLD, status, ier);
    else
      call MPI_Recv(func( first-1, t), 1, MPI_DOUBLE_PRECISION, id-1, id, MPI_COMM_WORLD, status, ier);
    end if
    

    if ( id==0 ) then
      call MPI_Send(func( first, t), 1, MPI_DOUBLE_PRECISION, np-1, id, MPI_COMM_WORLD, ier);
    else
      call MPI_Send(func( first, t), 1, MPI_DOUBLE_PRECISION, id-1, id, MPI_COMM_WORLD, ier);
    end if

    if ( id==np-1 ) then
      call MPI_Recv(func( last+1, t), 1, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, status, ier)
    else
      call MPI_Recv(func( last+1, t), 1, MPI_DOUBLE_PRECISION, id+1, id+1, MPI_COMM_WORLD, status, ier)
    end if

    

      do x=first,last
        func(x,t+1)=((func(x-1,t)+func(x+1,t))/2)-(G*(func(x+1,t)-func(x-1,t))/2)
      end do

  end do



  print *,11,id

  if ( id == 0 ) then
     
    param = 90
    open(14,file='Lax.dat')
    write(14,'(2(3x,e20.12))') (Real(x-1,8)*dx, func(x,param), x=first,last)
    close(14)

    call MPI_Send(param, 1, MPI_INTEGER, id+1, id+1, MPI_COMM_WORLD, ier)

    else
    call MPI_Recv(param, 1, MPI_INTEGER, id-1, id, MPI_COMM_WORLD, status, ier);

    open(14,file='Lax.dat', position="append")
!    do i = 1, 499
!      read(14,*)      
!    end do
    write(14,'(2(3x,e20.12))') (Real(x-1,8)*dx, func(x,param), x=first,last)
    close(14)

  end if
  call mpi_barrier(MPI_COMM_WORLD, ier)
  print *,"done"
  

  call MPI_FINALIZE(rc)
   
end
