module geometry

  implicit none

  interface angle
     real(8) function angle_xyz(x1, y1, z1, x2,y2,z2, x3,y3,z3)
       real(8) x1, y1, z1
       real(8) x2, y2, z2
       real(8) x3, y3, z3
     end function angle_xyz

     real(8) function angle_r(r1,r2,r3)
       real(8), dimension(1:3) :: r1,r2,r3
     endfunction angle_r
  end interface angle

  interface distance
     real(8) function distance_xyz(x1,y1,z1,x2,y2,z2)
       real(8) x1,y1,z1,x2,y2,z2
     end function distance_xyz

     real(8) function distance_r(r1,r2)
       real(8), dimension(1:3) :: r1,r2
     end function distance_r
  end interface distance
end module geometry


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(8) function angle_xyz(x1, y1, z1, x2,y2,z2, x3,y3,z3)
  real(8) x1, y1, z1
  real(8) x2, y2, z2
  real(8) x3, y3, z3

!		funct0ion, calculating angle A
!
!		A(x1,y1,z1) B(x2,y2,z2), C(x3,y3,z3)
!
!                      B
!                     /
!                    /
!                   /
!                  A -------------- C
!

  angle_xyz = acos( ( (x2-x1)*(x3-x1)+(y2-y1)*(y3-y1)+(z2-z1)*(z3-z1) )/ &
       (sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)*sqrt((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)) )
end function angle_xyz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(8) function distance_xyz(x1,y1,z1,x2,y2,z2)
  real(8) x1,y1,z1,x2,y2,z2
  distance_xyz = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
end function distance_xyz


real(8) function angle_r(r1,r2,r3)
  real(8), dimension(1:3) :: r1,r2,r3
  angle_r = acos((sum( (r2-r1)*(r3-r1))) / sqrt(sum((r2-r1)*(r2-r1))*sum((r3-r1)*(r3-r1))))
endfunction angle_r

real(8) function distance_r(r1,r2)
  real(8), dimension(1:3) :: r1, r2
  distance_r = sqrt(sum((r2-r1)*(r2-r1)))
endfunction distance_r
