module matrix_rotation

implicit none

	interface
	 integer function index_rotation_matrix(x,y,z)
		 real(8) x,y,z
	 end function index_rotation_matrix

	integer function cube2sphere(x, y, z, x_edge, y_edge, r_sphere, number_edge)
		real(8) x, y, z
		real(8) x_edge, y_edge, r_sphere
		integer number_edge
		real(8) rr
	end function cube2sphere
	end interface

CONTAINS



integer function init_compute_matrices(rots)
! computes matrices of rotation

	interface
		integer function cube2sphere(x,y,z,xe,ye,r,i)
			real(8) x,y,z,xe,ye,r
			integer i
		end function cube2sphere

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

		integer function index_rotation_matrix(x,y,z)
			real(8) x,y,z
		end function index_rotation_matrix
	end interface


	real(8), dimension(1:3,1:3, 48) :: rots

	real(8), parameter :: pi2 = 6283185307179586d-15
	real(8) rot(3,3), r(3)
	real(8) x,y,z, x_edge, y_edge
	integer i, j, index, status
	
	rots=0
	do i = 1,6
		do j = 1,12
			if (mod(j,3).ne.0) then
				x_edge = cos(pi2*dble(j)/12d0)
				y_edge = sin(pi2*dble(j)/12d0)
				status = cube2sphere(x, y, z, x_edge,y_edge, 1d0, i)

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




integer function index_rotation_matrix(x,y,z)
	interface
		 integer function index_verge_calc(x,y,z)
			 real(8) x,y,z
		 end function index_verge_calc

		 integer function index_vertex_calc(x,y,z)
			 real(8) x,y,z
		 end function index_vertex_calc

		 integer function matrix_verge_rotation(x,y,z,rot)
			 real(8) x,y,z, rot(1:3,1:3)
		 end function matrix_verge_rotation
	end interface

	real(8) x,y,z
	real(8) rot(1:3,1:3), r(1:3)
	integer status

	index_rotation_matrix = (index_verge_calc(x,y,z)-1) * 8
	status = matrix_verge_rotation(x,y,z,rot)
	r = matmul(rot, (/x,y,z/))
	index_rotation_matrix = index_rotation_matrix + index_vertex_calc(r(1),r(2),r(3))
	

end function index_rotation_matrix
