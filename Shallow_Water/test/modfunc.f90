	Module modfunc

	Use modnet

		Implicit None


		Public :: func

		Type func
			Integer(Kind=4) ns_x,nf_x,ns_y,nf_y

			Real(8), Allocatable :: d(:, :)
			Real(8), Allocatable :: du(:, :)
			Real(8), Allocatable :: dv(:, :)
		CONTAINS
			Procedure ::   init => func_init
			Procedure :: deinit => func_deinit
			Procedure :: eq => func_eq
		End Type

	CONTAINS

		Subroutine func_init(this,g)
			Class(func) :: this
			Class(grid) :: g

			call this.deinit()

			this.ns_x = g.ns_x-g.bstep
			this.nf_x = g.nf_x+g.fstep
			this.ns_y = g.ns_y-g.bstep
			this.nf_y = g.nf_y+g.bstep

			Allocate(this.d(this.ns_y:this.nf_y, this.ns_x:this.nf_x))
			Allocate(this.du(this.ns_y:this.nf_y, this.ns_x:this.nf_x))
			Allocate(this.dv(this.ns_y:this.nf_y, this.ns_x:this.nf_x))

		End Subroutine


		Subroutine func_deinit(this)
			Class(func) :: this

			if (Allocated(this.d)) Deallocate(this.d)
			if (Allocated(this.du)) Deallocate(this.du)
			if (Allocated(this.dv)) Deallocate(this.dv)

		End Subroutine


		Subroutine func_eq(this, that)

			Class(func) :: this, that
			Integer(4) x, y

			do y = this.ns_y, this.nf_y
				do x = this.ns_x,this.nf_x
					this.d(y,x)=that.d(y,x)
					this.du(y,x)=that.du(y,x)
					this.dv(y,x)=that.dv(y,x)
				end do
			end do

		End Subroutine


	End Module
