
integer function matrix_vertex_rotation(x,y,z, rot)
	use matmul_module
	use simple_rotations
	implicit none

	interface
		 integer function index_vertex_calc(x,y,z)
			 real(8) x,y,z
		 end function index_vertex_calc
	end interface

	real(8) x,y,z
	real(8) rot(1:3,1:3)

	integer index_vertex

	index_vertex = index_vertex_calc(x,y,z)
	matrix_vertex_rotation = 0

	select case (index_vertex)
		 case(1)
				rot=0d0; rot(1,1)=1d0; rot(2,2)=1d0; rot(3,3)=1d0!identify matrix
		 case(2)
				call matmul1(Rotation_xy, Rotation_mir, rot)
		 case(3)
				rot = Rotation_xy
		 case(4)
				call matmul1(Rotation_xy, Rotation_xy, Rotation_mir, rot)
		 case(5)
				call matmul1(Rotation_xy, Rotation_xy, rot)
		 case(6)
				call matmul1(transpose(Rotation_xy),Rotation_mir, rot)
		 case(7)
				rot = transpose(Rotation_xy)
		 case(8)
				rot = Rotation_mir
		 case default
				matrix_vertex_rotation = 1
	end select
	rot=transpose(rot)
	matrix_vertex_rotation = 1

end function matrix_vertex_rotation


integer function index_vertex_calc(x,y,z)
	real(8) x,y,z

	index_vertex_calc = 0

	If ( x .ge.  y  .and.    y  .ge.  0d0) then
		 index_vertex_calc = 1
	else if( y .ge.  x  .and.    x  .ge.  0d0) then
		 index_vertex_calc = 2
	else if( y .ge. -x  .and.   -x  .ge.  0d0) then
		 index_vertex_calc = 3
	else if(-x .ge.  y  .and.    y  .ge.  0d0) then
		 index_vertex_calc = 4
	else if(-x .ge. -y  .and.   -y  .ge.  0d0) then
		 index_vertex_calc = 5
	else if(-y .ge.  x  .and.   -x  .ge.  0d0) then
		 index_vertex_calc = 6
	else if(-y .ge.  x  .and.    x  .ge.  0d0) then
		 index_vertex_calc = 7
	else if( x .ge. -y  .and.   -y  .ge.  0d0) then
		 index_vertex_calc = 8
	end if

end function index_vertex_calc

integer function index_verge_calc(x,y,z)

	real(8) x,y,z

	If ( (-z) .ge. max(abs(x),abs(y)) ) then
		 index_verge_calc = 1
	else if(   x  .ge. max(abs(y),abs(z)) ) then
		 index_verge_calc = 2
	else if(   y  .ge. max(abs(x),abs(z)) ) then
		 index_verge_calc = 3
	else if(  -x  .ge. max(abs(y),abs(z)) ) then
		 index_verge_calc = 4
	else if(  -y  .ge. max(abs(x),abs(z)) ) then
		 index_verge_calc = 5
	else if(   z  .ge. max(abs(x),abs(y)) ) then
		 index_verge_calc = 6
	end if

end function index_verge_calc

