Module Schema

	Use modnet
	Use modfunc
	Use MPI
	Use method

IMPLICIT NONE

	Public :: sch

	Type sch

		Integer(4) casenumb, rcase, eqvtype

		Real(8), Allocatable :: k1(:)
		Real(8), Allocatable :: k2(:)
		Real(8), Allocatable :: k3(:)
		Real(8), Allocatable :: k4(:)


		CONTAINS

			Procedure, Public :: angle => sch_scheme_angle
			Procedure :: cross => sch_scheme_cross
			Procedure :: Lax => sch_scheme_Lax
			Procedure :: LaxVen => sch_scheme_LaxVen
			Procedure ::  RungeKutta=> sch_RungeKutta
			Procedure ::  FRunge=> sch_FRunge

			Procedure :: init => sch_init
			Procedure :: deinit => sch_deinit

	End Type

	CONTAINS

	Subroutine sch_init(this, casenumb, eqvtype, rcase, g)

		Class(grid) :: g
		Class(sch) :: this
		Integer(4), Intent(In) :: rcase, eqvtype, casenumb
		Integer(4) x

		this.casenumb = casenumb
		this.rcase = rcase
		this.eqvtype = eqvtype


		if ( this.casenumb == 5 ) then
			Allocate(this.k1(g.ns - g.bstep : g.nf + g.fstep))
			Allocate(this.k2(g.ns-g.bstep:g.nf+g.fstep))
			Allocate(this.k3(g.ns-g.bstep:g.nf+g.fstep))
			Allocate(this.k4(g.ns-g.bstep:g.nf+g.fstep))

			do x = g.ns, g.nf
				this.k1(x) = 0
				this.k2(x) = 0
				this.k3(x) = 0
				this.k4(x) = 0
			end do

! 			print *, "Allocate"
		end if

	End Subroutine

	Subroutine sch_scheme_angle(this, y, yprev, g)

		Class(sch) :: this
		Class(grid) :: g

		Integer(4) :: x
		Real(8), Intent(inout) :: y(g.ns - g.bstep : g.nf)
		Real(8), Intent(inout) :: yprev(g.ns - g.bstep : g.nf)
		Real(8) :: Gamma
		Gamma = g.Gamma

		do x=g.ns,g.nf
			if ( this.eqvtype == 2 ) then
				Gamma = (yprev(x)*g.dt)/g.dx
			end if
			y(x)=yprev(x)-Gamma*(yprev(x)-yprev(x-1))
			!y(x)=2
		end do


	End Subroutine


	Subroutine sch_scheme_cross(this, y, yprev, ynext, g)

		Class(sch) :: this
		Class(grid) :: g

		Integer(4) :: x
		Real(8), Intent(inout) :: y(g.ns - g.bstep : g.nf + g.fstep)
		Real(8), Intent(inout) :: yprev(g.ns - g.bstep : g.nf + g.fstep)
		Real(8), Intent(inout) :: ynext(g.ns - g.bstep : g.nf + g.fstep)
		Real(8) :: Gamma
		Gamma = g.Gamma

		do x=g.ns, g.nf
			if ( this.eqvtype == 2 ) then
				Gamma = (yprev(x)*g.dt)/g.dx
			end if
			ynext(x)=yprev(x)-Gamma*(y(x+1)-y(x-1))
		end do


	End Subroutine

	Subroutine sch_scheme_Lax(this, y, yprev, g)

		Class(sch) :: this
! 		Class(func) :: f, fprev
		Class(grid) :: g

		Integer(4) :: x
		Real(8), Intent(inout) :: y(g.ns - g.bstep : g.nf + g.fstep)
		Real(8), Intent(inout) :: yprev(g.ns - g.bstep : g.nf + g.fstep)
		Real(8) :: Gamma
		Gamma = g.Gamma

		do x=g.ns, g.nf
			if ( this.eqvtype == 2 ) then
				Gamma = (yprev(x)*g.dt)/g.dx
			end if
			y(x)=((yprev(x-1)+yprev(x+1))/2)-(Gamma*(yprev(x+1)-y(x-1))/2)
		end do


	End Subroutine


	Subroutine sch_scheme_LaxVen(this, f, fprev, fnext, g, m, t)

		Class(sch) :: this
		Class(func) :: f, fprev, fnext
		Class(grid) :: g
		Class(met) :: m

		Integer(4) :: x
		Integer(4) :: t

		if ( mod(t,2)==0  ) then
			call m.Message(f.d, g)
			do x=g.ns, g.nf
				call this.cross(f.d, fprev.d, fnext.d, g)
			end do
			call fprev.eq(fnext)
		else
			call m.Message(fprev.d, g)
			do x=g.ns, g.nf
					call this.Lax(f.d, fprev.d, g)
			end do
		end if

	End Subroutine



	Subroutine sch_RungeKutta(this, f, fprev, g, m)

		Class(sch) :: this
		Class(func) :: f, fprev
		Class(grid) :: g
		Class(met) :: m

		Integer(4) :: x
		Integer(4) :: ier
		Integer(4) :: status(MPI_STATUS_SIZE)
! 		Integer(4), Intent(In) :: eqvtype
		call MPI_Barrier(MPI_COMM_WORLD, ier)


		call this.FRunge(this.k1, fprev.d, g)
		call MPI_Barrier(MPI_COMM_WORLD, ier)
		call m.Message(this.k1, g)

		call this.FRunge(this.k2, fprev.d, g)
		call MPI_Barrier(MPI_COMM_WORLD, ier)
		call m.Message(this.k2, g)

		call this.FRunge(this.k3, fprev.d, g)
		call MPI_Barrier(MPI_COMM_WORLD, ier)
		call m.Message(this.k3, g)

		call this.FRunge(this.k4, fprev.d, g)
		call MPI_Barrier(MPI_COMM_WORLD, ier)
		call m.Message(this.k4, g)


		do x=g.ns, g.nf
			
			f.d(x) = fprev.d(x) + (this.k1(x) + 2.0*this.k2(x) + 2.0*this.k3(x) + this.k4(x))/6.0

		end do


	End Subroutine

	Subroutine sch_FRunge(this, k, y, g)
		Class(grid) :: g
		Class(sch) :: this

		Integer(4) :: x
! 		Integer(4), Intent(In) :: eqvtype

		Real(8), Intent(In) :: y(g.ns - g.bstep : g.nf + g.fstep)
		Real(8), Intent(inout) :: k(g.ns - g.bstep : g.nf + g.fstep)

			do x = g.ns, g.nf

				if ( this.rcase == 1 ) then
					k(x) = -g.velocity*(y(x)-y(x-1))
				elseif ( this.rcase == 2 ) then
					k(x) = -g.velocity*(y(x+1) - y(x-1))/2.0
				elseif ( this.rcase == 3 ) then
					k(x) = (- g.velocity*(8*(y(x+1) - y(x-1)) - (y(x+2) - y(x-2)))/12.0)
! 					k(x) = 1
				end if

			end do


	end Subroutine


	Subroutine sch_deinit(this)
		Class(sch) :: this

		if (Allocated(this.k1)) Deallocate(this.k1)
		if (Allocated(this.k2)) Deallocate(this.k2)
		if (Allocated(this.k3)) Deallocate(this.k3)
		if (Allocated(this.k4)) Deallocate(this.k4)

	End Subroutine


End Module