program uravnenie

  Use MPI

  Implicit None


  real(8) lengthX, lengthT, u, dx, dt
  integer x,t, ier, id,np, parallparam, parallparamplus, synchr, i, rc, parnumb
  integer Xmax, Tmax, param
  real(8), Allocatable :: func(:,:)

  integer status(MPI_STATUS_SIZE)  
  

  call MPI_Init(ier)

  call MPI_Comm_rank(MPI_COMM_WORLD,id,ier)
  call MPI_Comm_size(MPI_COMM_WORLD,np,ier)
  

  lengthX=10
  lengthT=10
  Xmax = 1000
  Tmax = 1000
  dx = lengthX/Xmax; dt = lengthT/Tmax
  u=.1
  parallparam = id*Xmax/np
  parallparamplus = (id+1)*Xmax/np 
  synchr = 0
  parnumb = Xmax - parallparamplus


  write(*,'(2x,a,i4,3x,a,i4,3x,a,f14.7)') "myid =",id,"np=",np,"dx=",dx

  Allocate(func(Xmax, Tmax))
  
   if ( id == 0) then
      open(14,file ='func.dat')
      do x = 1,Xmax
        func(x,1) = exp(-0.0002*(x-500)**2)
        write(14,'(2(3x,e20.12))') Real(x-1,8)*dx, func(x,1)
      end do
      close(14)


      ! write(*,*)"Choose the type of differenceses you want to use: ", &
      ! "21 - forvard t backward x, 22 - forvard t forvard x, 23 - forvard t central x", &
      ! "31 - central t backward x, 32 - central t forvard x, 33 - central t central x"

    end if

  param = 21
   
  do t=1,Tmax
    func(1,t)=func(Xmax, t)
  end do

  do x=2,Xmax
    func(x,2)=func(x,1)-u*(dt/dx)*(func(x,1)-func((x-1),1))
  end do
     

  if ( param == 21) then

    if ( id == 0)then

      synchr=1
      do i=1,np-1
        call MPI_Send(synchr, 1, MPI_INTEGER,i, 1, MPI_COMM_WORLD, ier);
      end do

    else

      call MPI_Recv(synchr,1,MPI_INTEGER, 0, 1, MPI_COMM_WORLD, status, ier);
      write(*, *) "synchr=" , synchr     
      
    end if 



      do t=1,Tmax-1

        if ( id == 0 ) then
          do x=2,parallparamplus
              func(x,(t+1))=func(x,t)-u*(dt/dx)*(func(x,t)-func((x-1),t))
              call MPI_Recv(func(parallparamplus, t), parnumb, MPI_REAL, 1, 0, MPI_COMM_WORLD, status, ier);
          end do

        else
          do x=parallparam,Xmax
            func(x,(t+1))=func(x,t)-u*(dt/dx)*(func(x,t)-func((x-1),t))
            call MPI_Send(func(parallparam, t), parnumb, MPI_REAL, 0, 0, MPI_COMM_WORLD, ier);
          end do
        end if

      end do

  elseif ( param == 22 ) then

      do t=1,Tmax-1
        do x=2,Xmax
          if(x<Xmax)then
            func(x,(t+1))=func(x,t)-u*(dt/dx)*(func((x+1),t)-func(x,t))
          else
            func(x,(t+1))=func(x,t)-u*(dt/dx)*(func(1,t)-func(x,t))
          end if
        end do
      end do

  elseif ( param == 23 ) then

      do t=1,Tmax-1
        do x=2,Xmax
          if(x<Xmax)then
            func(x,(t+1))=func(x,t)-u*(dt/dx/2)*(func((x+1),t)-func(x-1,t))
          else
            func(x,(t+1))=func(x,t)-u*(dt/dx/2)*(func(2,t)-func(x-1,t))
          end if
        end do
      end do

  elseif ( param == 31 ) then

      do t=2,Tmax-1
        do x=2,Xmax
            func(x,(t+1))=func(x,t-1)-u*(2*dt/dx)*(func(x,t)-func((x-1),t))
        end do
      end do

  elseif ( param == 32 ) then

      do t=2,Tmax-1
        do x=2,Xmax
          if(x<Xmax)then
            func(x,(t+1))=func(x,t-1)-u*(2*dt/dx)*(func((x+1),t)-func(x,t))
          else
            func(x,(t+1))=func(x,t-1)-u*(2*dt/dx)*(func(1,t)-func(x,t))
          end if
        end do
      end do

  elseif ( param == 33 ) then

      do t=2,Tmax-1
        do x=2,Xmax
          if(x<Xmax)then
            func(x,(t+1))=func(x,t-1)-u*(dt/dx)*(func((x+1),t)-func(x-1,t))
          else
            func(x,(t+1))=func(x,t-1)-u*(dt/dx)*(func(2,t)-func(x-1,t))
          end if
        end do
      end do

  end if
     



    if ( id == 0 ) then
      !write(*,*) "Choose the moment you want to see: " 
      param = 10
      
      open(14,file='func1.dat')
      
      write(14,'(2(3x,e20.12))') (Real(x-1,8)*dx, func(x,param), x=1,Xmax)
      
      close(14)
    end if

    if ( id == 1 ) then
      !write(*,*) "Choose the moment you want to see: " 
      param = 10
      
      open(14,file='func2.dat')
      
      write(14,'(2(3x,e20.12))') (Real(x-1,8)*dx, func(x,param), x=1,Xmax)
      
      close(14)
    end if

         ! z = Real(i,8)
         ! floor mod round

  call MPI_FINALIZE(rc)
   
end
