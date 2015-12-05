program uravnenie

  Implicit None

  real(8) lengthX, lengthT, u, dx, dt
  integer x,t
  integer Xmax, Tmax
  real, Allocatable :: func(:,:)  
  
  lengthX=10
  lengthT=10
  Xmax = 1000
  Tmax = 1000
  dx = lengthX/Xmax; dt = lengthT/Tmax
  u=.2
  
  Allocate(func(Xmax, Tmax))
  
   
    open(14,file='func.dat')
   do x = 1,Xmax
    func(x,1) = exp(-0.0002*(x-500)**2)
    
    write(14,'(2(3x,e20.12))') Real(x-1,8)*dx, func(x,1)
!     print *,x,func(x,1)
   end do
   close(14)
   
     do t=1,Tmax
     func(1,t)=func(Xmax, t)
     end do
     
     do t=1,Tmax-1
        if(t==95)then
          open(14,file='func1.dat')
          do x=2,Xmax
            if(x<Xmax)then
              func(x,(t+1))=func(x,t)-u*(dt/dx)*(func((x+1),t)-func(x,t))
            else
              func(x,(t+1))=func(x,t)-u*(dt/dx)*(func(1,t)-func(x,t))
            end if
            write(14,'(2(3x,e20.12))') Real(x-1,8)*dx, func(x,t)
          end do
          close(14)
        else
          do x=2,Xmax
            if(x<Xmax)then
              func(x,(t+1))=func(x,t)-u*(dt/dx)*(func((x+1),t)-func(x,t))
            else
              func(x,(t+1))=func(x,t)-u*(dt/dx)*(func(1,t)-func(x,t))
            end if
          end do
        end if
     end do
    
    STOP
   
   end
  
!  contains
!  real function f(x, t)
!  integer x, t
  

 !   f(x,t)=f(x,(t-1)) - u*(dt/dx)*(f(x+1,t-1) - f(x,t-1))
   