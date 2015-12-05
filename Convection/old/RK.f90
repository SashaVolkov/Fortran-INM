Module RungeKutta

Use modnet, modfunc, MPI


IMPLICIT NONE

	Public :: rk

	Type met

		Real(8) k1, k2, k3, k4


End Module