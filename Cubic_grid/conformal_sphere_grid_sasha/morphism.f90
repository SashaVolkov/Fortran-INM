module morphism

interface
	 integer function cube2sphere(x, y, z, x_edge, y_edge, r_sphere, number_edge)
		 real(8) x, y, z
		 real(8) x_edge, y_edge, r_sphere
		 integer number_edge
	 end function cube2sphere

	 integer function reverse_stereo_projection(x_plane, y_plane, r_sphere, x, y, z)
		 real(8) x_plane, y_plane, r_sphere, x, y, z
	 end function reverse_stereo_projection

	 integer function stereo_projection(x_plane, y_plane, r_sphere, x, y, z)
		 real(8) x_plane, y_plane, r_sphere, x, y, z
	 end function stereo_projection
	 
	 integer function matrix_rotation_to_top(x,y,z, rot)
		 real(8) x,y,z
		 real(8), dimension(1:3,1:3) :: rot
	 end function matrix_rotation_to_top

	 integer function matrix_verge_rotation(x,y,z, rot)
		 real(8) x,y,z
		 real(8), dimension(1:3,1:3) :: rot
	 end function matrix_verge_rotation

end interface
end module morphism

!functions2