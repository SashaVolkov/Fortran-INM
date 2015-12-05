	Module modfunc


		Implicit None


		Public :: func

		Type func
			Integer(Kind=4) ns,nf

			Real(8), Allocatable :: d(:)
		CONTAINS
			Procedure ::   init => func_init
			Procedure :: deinit => func_deinit
			Procedure :: eq => func_eq
		End Type

	CONTAINS

		Subroutine func_init(this,ns,nf)
			Class(func) :: this
			Integer(Kind=4), Intent(In) :: ns,nf

			call this.deinit()

			this.ns = ns
			this.nf = nf

			Allocate(this.d(this.ns:this.nf))

		End Subroutine




		Subroutine func_deinit(this)
			Class(func) :: this

			if (Allocated(this.d)) Deallocate(this.d)

		End Subroutine



		Subroutine func_eq(this, that)

			Class(func) :: this, that
			Integer(4) i

			do i = this.ns,this.nf
				this.d(i)=that.d(i)  	
			end do

		End Subroutine





	End Module


