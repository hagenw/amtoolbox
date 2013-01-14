SUBROUTINE InitializeMiddleEar
	USE Declare
	IMPLICIT NONE

	REAL(8) stapes_area
	stapes_area = stapesArea(parameterSet)

	   !the ME represents a resistance Rme
           !that equals the cochlear input impedance at low frequencies. 
           !By doing this, there are little reflections of cochlear energy
           !back into the model, as it absorbs the energy at the stapes.
                 
                  !Rme=SQRT(ZweigMso * Ko) !kg/m4s
                  Mme=1d0 !test for now !!was 1
                  !debug(1:n)= Rme
                  !debug(2)= Ko
                  !debug(3)= Mme
                  !d_m_factor is going to be multiplied with Vbm0 for g(0) in RK4
	        q0_factor = ZweigMpo * bm_width
                p0x = (ZweigMso * dx) / (Mme * ZweigMpo * bm_width )
                d_m_factor = - p0x * stapes_area * Rme
                !debug(4)=d_m_factor
                !d_m_factor = ( Rme * dx ) / ZweigMpo
                RK4_0 = -(bm_width * ZweigMpo)/(Mme * stapes_area) !for RK4 method 
                RK4G_0= (ZweigMpo * bm_width)/(ZweigMso * stapes_area * dx)
        !All the other cases, old implementation and unchanged (sv)
END SUBROUTINE InitializeMiddleEar
