SUBROUTINE InitializeZweig

	USE Declare
	IMPLICIT NONE
	REAL(8) f_resonance(1:n)
	REAL(8) exact_delay(1:n)
	REAL(8) relative_position(n)
	REAL(8) adaptation(n)
	REAL(8) HighFrequencyAdaptationFactor
               

        !Delay corresponds to the maximal possible delay in the model
        !determines the length of the buffers
        f_resonance = omega(1:n) / (2d0 * pi)
        exact_delay = SheraMuMax / (f_resonance(1:n) * dt) ![]
        !exact_delay = ZweigFactor / (f_resonance(1:n) * dt) ![]
        ! exact delay in (fraction of) number of samples
        delay = INT (exact_delay) + 1.d0 ! delay in integer number of samples
           !debug(1:n)=0
           !debug(1:n)=delay(1:n)
           !delay_deviation = DFLOAT(delay) - exact_delay
           !!comment this out when going with exact buffer!
 
        ALLOCATE (Ybuffer(SUM(delay)), STAT = err)
        IF (err/=0) CALL AllocationError

 	Ybuffer = 0
        ZweigSample1(1) = 1d0 !write pointer
        Zwp(1)=1d0            !write pointer
        DO i = 2, n
           ZweigSample1(i) = ZweigSample1(i-1) + delay(i-1)
           Zwp(i)=Zwp(i-1) + delay(i-1)
        ENDDO
           ZweigSample2 = ZweigSample1 + 1d0    

END SUBROUTINE InitializeZweig
