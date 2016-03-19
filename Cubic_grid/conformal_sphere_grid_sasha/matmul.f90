module matmul_module
implicit none
interface matmul1
   module procedure matmul2
   module procedure matmul3
   module procedure matmul4
   module procedure matmul5
   module procedure matmul6
end interface matmul1
CONTAINS

	subroutine matmul2(A,B,Z)
		real*8 A(:,:),B(:,:),Z(:,:)
		Z=matmul(A,B)
	end subroutine

	subroutine matmul3(A,B,C,Z)
		real*8 A(:,:),B(:,:),C(:,:), Z(:,:)
		Z = matmul(A,matmul(B,C))
	end subroutine

	subroutine matmul4(A,B,C,D,Z)
		real*8 A(:,:),B(:,:),C(:,:),D(:,:), Z(:,:)
		Z = matmul(A, matmul(B, matmul(C,D)))
	end subroutine


	subroutine matmul5(A,B,C,D,E, Z)
		real*8 A(:,:),B(:,:),C(:,:),D(:,:),E(:,:),Z(:,:)
		Z = matmul(A, matmul(B, matmul(C,matmul(D,E))))
	end subroutine


	subroutine matmul6(A,B,C,D,E,F,Z)
		real*8 A(:,:),B(:,:),C(:,:),D(:,:),E(:,:),F(:,:), Z(:,:)
		Z = matmul(A, matmul(B, matmul(C,matmul(D,matmul(E,F)))))
	end subroutine
end module matmul_module
