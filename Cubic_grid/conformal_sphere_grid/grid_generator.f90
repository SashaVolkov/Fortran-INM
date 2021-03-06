module grid_generator
	interface 
		 integer function conformal_cubed_sphere_grid_generation(x_points, y_points)
			 integer x_points, y_points
		 end function conformal_cubed_sphere_grid_generation
		 integer function cubed_sphere_grid_generation(x_points, y_points)
			 integer x_points, y_points
		 end function cubed_sphere_grid_generation
	endinterface
endmodule

integer function conformal_cubed_sphere_grid_generation(x_points,y_points)
	use morphism
	use matrix_rotation
	use conformal

	implicit none
	integer x_points, y_points
	character*14 filename
	character istring

	real*8 rots(3,3,48)
	real*8 rot(1:3,1:3)
	real*8 r(1:3)
	real*8 x,y,z
	real*8 x1,y1, x_edge, y_edge
	complex*16 w

	integer i, j, k, status, index
	status = init_compute_matrices(rots)
	
	open (20, file = "grid/matrices.dat")
	do i =1,48
		 write(20,*) rots(1:3,1,i)
		 write(20,*) rots(1:3,2,i)
		 write(20,*) rots(1:3,3,i)
		 write(20,*)
	end do
	close(20)

	do i=1, 6
		write(istring(1:1), '(i1.1)') i
		filename = "grid/edge" // istring // ".dat" 
		open (20, file = filename)
		do j=-x_points,x_points
			do k=-y_points,y_points
				x1=j/dble(x_points)
				y1=k/dble(y_points)

				status = cube2sphere(x,y,z,x1,y1,1d0,i)
				index = index_rotation_matrix(x,y,z)

				if (abs(x1)>abs(y1)) then
					y_edge = abs(sign(1d0,x1)*(1-abs(x1)))/2d0
					x_edge = abs(sign(1d0,y1)*(1-abs(y1)))/2d0
				else
					x_edge = abs(sign(1d0,x1)*(1-abs(x1)))/2d0
					y_edge = abs(sign(1d0,y1)*(1-abs(y1)))/2d0
				endif

				status = conformal_z_w(dcmplx(x_edge,y_edge), w)

				status = reverse_stereo_projection(dreal(w),dimag(w),1d0,x,y,z)

				r = matmul(transpose(rots(1:3,1:3,index)),(/x,y,z/))
				x = r(1); y=r(2); z=r(3)


				write(20,*) x,y,z

			end do
			write(20,*)
		end do
	end do
	close(20)

	open(20, file = "grid/parameters.dat")
	write(20,*) x_points, y_points
	close(20)
end function conformal_cubed_sphere_grid_generation

integer function cubed_sphere_grid_generation(x_points, y_points)
	use morphism
	use matrix_rotation

	implicit none
	integer x_points, y_points
	character*14 filename
	character istring

	real*8 x,y,z
	real*8 x1,y1

	integer i, j, k, status
	do i=1, 6
		 write(istring(1:1), '(i1.1)') i
		 filename = "grid/edge" // istring // ".dat" 
		 open (20, file = filename)
		 do j=-x_points,x_points
				do k=-y_points,y_points
					 x1=j/dble(x_points)
					 y1=k/dble(y_points)
					 
					 status = cube2sphere(x,y,z,x1,y1,1d0,i)
					 write(20,*) x,y,z
					 
				end do
		 end do
	end do
	close(20)

	open(20, file = "grid/parameters.dat")
	write(20,*) x_points, y_points
	close(20)


end function cubed_sphere_grid_generation
