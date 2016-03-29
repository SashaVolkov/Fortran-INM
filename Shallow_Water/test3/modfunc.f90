	Module modfunc

	Use modnet, Only: grid

		Implicit None

		Private
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

			this.ns_x = g.first_x
			this.nf_x = g.last_x
			this.ns_y = g.first_y
			this.nf_y = g.last_y

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
