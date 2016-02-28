Module Schema

	Use modnet
	Use MPI
	Use method

IMPLICIT NONE

	Public :: sch

	Type sch

		Integer(4) casenumb, rcase, eqvtype
		Real(8) grav, cor, height

		Real(8), Allocatable :: ku1(:,:)
		Real(8), Allocatable :: ku2(:,:)
		Real(8), Allocatable :: ku3(:,:)
		Real(8), Allocatable :: ku4(:,:)

		Real(8), Allocatable :: kv1(:,:)
		Real(8), Allocatable :: kv2(:,:)
		Real(8), Allocatable :: kv3(:,:)
		Real(8), Allocatable :: kv4(:,:)

		Real(8), Allocatable :: kh1(:,:)
		Real(8), Allocatable :: kh2(:,:)
		Real(8), Allocatable :: kh3(:,:)
		Real(8), Allocatable :: kh4(:,:)

		CONTAINS

			Procedure, Public :: linear => sch_scheme_linear
! 			Procedure :: cross => sch_scheme_cross
! 			Procedure :: Lax => sch_scheme_Lax
! 			Procedure :: LaxVen => sch_scheme_LaxVen
			Procedure ::  RungeKutta=> sch_RungeKutta
			Procedure ::  FRunge=> sch_FRunge

			Procedure :: init => sch_init
			Procedure :: deinit => sch_deinit

	End Type

	CONTAINS

	Subroutine sch_init(this, casenumb, eqvtype, rcase, cor, grav, height, g)

		Class(grid) :: g
		Class(sch) :: this
		Integer(4), Intent(In) :: rcase, eqvtype, casenumb
		Real(8), Intent(In) ::  cor, grav, height
		Integer(4) x

		this.casenumb = casenumb
		this.rcase = rcase
		this.eqvtype = eqvtype
		this.grav = grav
		this.cor = cor
		this.height = height


		if ( this.casenumb == 5 ) then
			Allocate(this.ku1(g.ns_y - g.bstep : g.nf_y + g.fstep, g.ns_x - g.bstep : g.nf_x + g.fstep))
			Allocate(this.ku2(g.ns_y - g.bstep : g.nf_y + g.fstep, g.ns_x - g.bstep : g.nf_x + g.fstep))
			Allocate(this.ku3(g.ns_y - g.bstep : g.nf_y + g.fstep, g.ns_x - g.bstep : g.nf_x + g.fstep))
			Allocate(this.ku4(g.ns_y - g.bstep : g.nf_y + g.fstep, g.ns_x - g.bstep : g.nf_x + g.fstep))

			Allocate(this.kv1(g.ns_y - g.bstep : g.nf_y + g.fstep, g.ns_x - g.bstep : g.nf_x + g.fstep))
			Allocate(this.kv2(g.ns_y - g.bstep : g.nf_y + g.fstep, g.ns_x - g.bstep : g.nf_x + g.fstep))
			Allocate(this.kv3(g.ns_y - g.bstep : g.nf_y + g.fstep, g.ns_x - g.bstep : g.nf_x + g.fstep))
			Allocate(this.kv4(g.ns_y - g.bstep : g.nf_y + g.fstep, g.ns_x - g.bstep : g.nf_x + g.fstep))

			Allocate(this.kh1(g.ns_y - g.bstep : g.nf_y + g.fstep, g.ns_x - g.bstep : g.nf_x + g.fstep))
			Allocate(this.kh2(g.ns_y - g.bstep : g.nf_y + g.fstep, g.ns_x - g.bstep : g.nf_x + g.fstep))
			Allocate(this.kh3(g.ns_y - g.bstep : g.nf_y + g.fstep, g.ns_x - g.bstep : g.nf_x + g.fstep))
			Allocate(this.kh4(g.ns_y - g.bstep : g.nf_y + g.fstep, g.ns_x - g.bstep : g.nf_x + g.fstep))
		end if

	End Subroutine

	Subroutine sch_scheme_linear(this, h, hpr, u, upr, v, vpr, g)

		Class(sch) :: this
		Class(grid) :: g

		Integer(4) :: x, y
		Real(8), Intent(out) :: h(g.ns_y - g.bstep : g.nf_y + g.fstep, g.ns_x - g.bstep : g.nf_x + g.fstep)
		Real(8), Intent(in) :: hpr(g.ns_y - g.bstep : g.nf_y + g.fstep, g.ns_x - g.bstep : g.nf_x + g.fstep)
		Real(8), Intent(out) :: u(g.ns_y - g.bstep : g.nf_y + g.fstep, g.ns_x - g.bstep : g.nf_x + g.fstep)
		Real(8), Intent(in) :: upr(g.ns_y - g.bstep : g.nf_y + g.fstep, g.ns_x - g.bstep : g.nf_x + g.fstep)
		Real(8), Intent(out) :: v(g.ns_y - g.bstep : g.nf_y + g.fstep, g.ns_x - g.bstep : g.nf_x + g.fstep)
		Real(8), Intent(in) :: vpr(g.ns_y - g.bstep : g.nf_y + g.fstep, g.ns_x - g.bstep : g.nf_x + g.fstep)
		Real(8) cor, gr, height

		cor = this.cor; gr = this.grav; height = this.height

		do y = g.ns_y, g.nf_y
			do x = g.ns_x, g.nf_y
				u(y,x) = upr(y,x) + g.dt*(cor*vpr(y,x) - (gr*0.5)*(hpr(y, x+1) - hpr(y, x-1))/g.dx)! - (uprev(y,x)*(uprev(y, x+1) - uprev(y, x-1)) + vprev(y,x)*(uprev(y+1, x) - uprev(y-1, x)))
				v(y,x) = vpr(y,x) - g.dt*(cor*upr(y,x) - (gr*0.5)*(hpr(y+1, x) - hpr(y-1, x))/g.dx)! - (uprev(y,x)*(vprev(y, x+1) - vprev(y, x-1)) + vprev(y,x)*(vprev(y+1, x) - vprev(y-1, x)))
				h(y,x) = hpr(y,x) - g.dt*height*0.5*(((upr(y, x+1) - upr(y, x-1)))/g.dx + ((vpr(y+1, x) - vpr(y-1, x)))/g.dy)! - (uprev(y,x)*(hprev(y, x+1) - hprev(y, x-1)) + vprev(y,x)*(hprev(y+1, x) - hprev(y-1, x)))
			end do
		end do
	End Subroutine



	Subroutine sch_RungeKutta(this, f, fprev, g, m)

		Class(sch) :: this
		Class(func) :: f, fprev
		Class(grid) :: g
		Class(met) :: m

		Integer(4) :: x,y
		Integer(4) :: ier
		Integer(4) :: status(MPI_STATUS_SIZE)
! 		Integer(4), Intent(In) :: eqvtype
		call MPI_Barrier(MPI_COMM_WORLD, ier)


		call this.FRunge(this.ku1, this.kv1, this.kh1, fprev.du, fprev.dv, fprev.d, g)
		call MPI_Barrier(MPI_COMM_WORLD, ier)
! 		call m.Message(this.ku1, g);call m.Message(this.kv1, g);call m.Message(this.kh1, g)

		call this.FRunge(this.ku2, this.kv2, this.kh2, fprev.du, fprev.dv, fprev.d, g)
		call MPI_Barrier(MPI_COMM_WORLD, ier)
! 		call m.Message(this.ku2, g);call m.Message(this.kv2, g);call m.Message(this.kh2, g)

		call this.FRunge(this.ku3, this.kv3, this.kh3, fprev.du, fprev.dv, fprev.d, g)
		call MPI_Barrier(MPI_COMM_WORLD, ier)
! 		call m.Message(this.ku3, g);call m.Message(this.kv3, g);call m.Message(this.kh3, g)

		call this.FRunge(this.ku4, this.kv4, this.kh4, fprev.du, fprev.dv, fprev.d, g)
		call MPI_Barrier(MPI_COMM_WORLD, ier)
! 		call m.Message(this.ku4, g);call m.Message(this.kv4, g);call m.Message(this.kh4, g)

		do y = g.ns_y, g.nf_y
			do x=g.ns_x, g.nf_x
				f.du(y,x) = fprev.du(y,x) + (this.ku1(y,x) + 2.0*this.ku2(y,x) + 2.0*this.ku3(y,x) + this.ku4(y,x))/6.0
				f.dv(y,x) = fprev.dv(y,x) + (this.kv1(y,x) + 2.0*this.kv2(y,x) + 2.0*this.kv3(y,x) + this.kv4(y,x))/6.0
				f.d(y,x) = fprev.d(y,x) + (this.kh1(y,x) + 2.0*this.kh2(y,x) + 2.0*this.kh3(y,x) + this.kh4(y,x))/6.0
			end do
		end do

	End Subroutine


	Subroutine sch_FRunge(this, ku, kv, kh, upr, vpr, hpr, g)
		Class(grid) :: g
		Class(sch) :: this

		Integer(4) :: x,y
! 		Integer(4), Intent(In) :: eqvtype

		Real(8), Intent(In) :: upr(g.first_y : g.last_y, g.first_x : g.last_x)
		Real(8), Intent(In) :: vpr(g.first_y : g.last_y, g.first_x : g.last_x)
		Real(8), Intent(In) :: hpr(g.first_y : g.last_y, g.first_x : g.last_x)
		Real(8), Intent(inout) :: ku(g.first_y : g.last_y, g.first_x : g.last_x)
		Real(8), Intent(inout) :: kv(g.first_y : g.last_y, g.first_x : g.last_x)
		Real(8), Intent(inout) :: kh(g.first_y : g.last_y, g.first_x : g.last_x)

		Real(8) cor, gr, height
		cor = this.cor; gr = this.grav; height = this.height

		do y = g.ns_y, g.nf_y
			do x = g.ns_x, g.nf_x
				ku(y,x) = g.dt*(cor*vpr(y,x) - (gr*0.5)* &
					(4.0*hpr(y, x+1) + 0.5*hpr(y, x-2) - 4.0*hpr(y, x-1) - 0.5*hpr(y, x+2))/(g.dx*3.0))


				kv(y,x) = g.dt*( - cor*upr(y,x) - (gr*0.5)* &
					(4.0*hpr(y+1, x) + 0.5*hpr(y-2, x) - 4.0*hpr(y-1, x) - 0.5*hpr(y+2, x))/(g.dy*3.0))


				kh(y,x) = - g.dt*height*0.5*(&
					4.0*((upr(y, x+1) - upr(y, x-1))/(g.dx*3.0) + (vpr(y+1, x) - vpr(y-1, x))/(g.dy*3.0)) - &
					0.5*((upr(y, x+2) - upr(y, x-2))/(g.dx*3.0) + (vpr(y+2, x) - vpr(y-2, x))/(g.dy*3.0)))

			end do
		end do

	end Subroutine





	Subroutine sch_deinit(this)
		Class(sch) :: this

		if (Allocated(this.ku1)) Deallocate(this.ku1)
		if (Allocated(this.ku2)) Deallocate(this.ku2)
		if (Allocated(this.ku3)) Deallocate(this.ku3)
		if (Allocated(this.ku4)) Deallocate(this.ku4)

		if (Allocated(this.kv1)) Deallocate(this.kv1)
		if (Allocated(this.kv2)) Deallocate(this.kv2)
		if (Allocated(this.kv3)) Deallocate(this.kv3)
		if (Allocated(this.kv4)) Deallocate(this.kv4)

		if (Allocated(this.kh1)) Deallocate(this.kh1)
		if (Allocated(this.kh2)) Deallocate(this.kh2)
		if (Allocated(this.kh3)) Deallocate(this.kh3)
		if (Allocated(this.kh4)) Deallocate(this.kh4)

	End Subroutine


End Module






! 	Subroutine sch_scheme_cross(this, f, fprev, fnext, g)

! 		Class(sch) :: this
! 		Class(grid) :: g

! 		Integer(4) :: x
! 		Real(8), Intent(inout) :: f(g.ns - g.bstep : g.nf + g.fstep)
! 		Real(8), Intent(inout) :: fprev(g.ns - g.bstep : g.nf + g.fstep)
! 		Real(8), Intent(inout) :: fnext(g.ns - g.bstep : g.nf + g.fstep)
! 		Real(8) :: Gamma
! 		Gamma = g.Gamma

! 		do x=g.ns, g.nf
! 			if ( this.eqvtype == 2 ) then
! 				Gamma = (fprev(x)*g.dt)/g.dx
! 			end if
! 			fnext(x)=fprev(x)-Gamma*(f(x+1)-f(x-1))
! 		end do


! 	End Subroutine

! 	Subroutine sch_scheme_Lax(this, f, fprev, g)

! 		Class(sch) :: this
! 		Class(grid) :: g

! 		Integer(4) :: x
! 		Real(8), Intent(inout) :: f(g.ns - g.bstep : g.nf + g.fstep)
! 		Real(8), Intent(inout) :: fprev(g.ns - g.bstep : g.nf + g.fstep)
! 		Real(8) :: Gamma
! 		Gamma = g.Gamma

! 		do x=g.ns, g.nf
! 			if ( this.eqvtype == 2 ) then
! 				Gamma = (fprev(x)*g.dt)/g.dx
! 			end if
! 			f(x)=((fprev(x-1)+fprev(x+1))/2)-(Gamma*(fprev(x+1)-f(x-1))/2)
! 		end do


! 	End Subroutine


! 	Subroutine sch_scheme_LaxVen(this, f, fprev, fnext, g, m, t)

! 		Class(sch) :: this
! 		Class(func) :: f, fprev, fnext
! 		Class(grid) :: g
! 		Class(met) :: m

! 		Integer(4) :: x
! 		Integer(4) :: t

! 		if ( mod(t,2)==0  ) then
! 			call m.Message(f.d, g)
! 			do x=g.ns, g.nf
! 				call this.cross(f.d, fprev.d, fnext.d, g)
! 			end do
! 			call fprev.eq(fnext)
! 		else
! 			call m.Message(fprev.d, g)
! 			do x=g.ns, g.nf
! 					call this.Lax(f.d, fprev.d, g)
! 			end do
! 		end if

! 	End Subroutine
