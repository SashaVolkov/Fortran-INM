module matrix_rotation

use subfunc, Only: func

implicit none

	interface
	 integer function index_rotation_matrix(x,y,z)
		 real(8) x,y,z
	 end function index_rotation_matrix
	end interface

CONTAINS



	integer function init_compute_matrices(rots)
	! computes matrices of rotation

		interface
			integer function index_rotation_matrix(x,y,z)
				real(8) x,y,z
			end function index_rotation_matrix

			integer function matrix_verge_rotation(x,y,z,rot)
				real(8) x,y,z
				real(8) rot(3,3)
			end function matrix_verge_rotation

			integer function matrix_vertex_rotation(x,y,z,rot)
				real(8) x,y,z
				real(8) rot(3,3)
			end function matrix_vertex_rotation

			integer function matrix_rotation_to_top(x,y,z,rot)
				real(8) x,y,z
				real(8) rot(3,3)
			end function matrix_rotation_to_top
		end interface


		real(8), dimension(1:3,1:3, 48) :: rots

		real(8), parameter :: pi2 = 6283185307179586d-15
		real(8) rot(3,3), r(3)
		real(8) x,y,z, x_edge, y_edge
		integer i, j, index, status
		Type(func) :: f
		
		rots=0
		do i = 1,6
			do j = 1,12
				if (mod(j,3).ne.0) then
					x_edge = cos(pi2*dble(j)/12d0)
					y_edge = sin(pi2*dble(j)/12d0)
					call f.cube2sphere(x, y, z, x_edge,y_edge, 1d0, i, status)

					index = index_rotation_matrix(x,y,z)

					status = matrix_verge_rotation(x,y,z, rot)
					rots(1:3,1:3,index) = rot

					r = matmul(rot,(/x,y,z/))
					x=r(1);y=r(2);z=r(3)

					if ( matrix_vertex_rotation(x,y,z,rot) ) then
						rots(1:3,1:3,index) = matmul(rot, rots(1:3,1:3,index))
					end if 

					r = matmul(rot,(/x,y,z/))
					x=r(1); y=r(2); z=r(3)
					 
					if ( matrix_rotation_to_top(x,y,z,rot) ) then
						rots(1:3,1:3,index) = matmul(rot, rots(1:3,1:3,index))
					end if 

				end if
			end do
		end do
		init_compute_matrices = 1
	end function init_compute_matrices


end module matrix_rotation


