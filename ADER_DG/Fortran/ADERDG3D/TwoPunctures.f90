!
! This is a Fortran adapter to the TwoPunctures C++ code.
! It makes use of the plain C TwoPuncturesAdapter.c.
!
! This is only needed when you want to use TwoPunctures from your
! Fortran code.
!
#ifdef TWOPUNCTURES 
MODULE TwoPunctures_mod 
	IMPLICIT NONE

	INTERFACE
		SUBROUTINE TwoPunctures_Setup() BIND(C)
			USE, INTRINSIC :: ISO_C_BINDING
			IMPLICIT NONE
		END SUBROUTINE TwoPunctures_Setup
	END INTERFACE
	
	INTERFACE
		SUBROUTINE TwoPunctures_AmendParameters() BIND(C)
			USE, INTRINSIC :: ISO_C_BINDING
			IMPLICIT NONE
		END SUBROUTINE TwoPunctures_AmendParameters
	END INTERFACE

	INTERFACE
		SUBROUTINE TwoPunctures_Run() BIND(C)
			USE, INTRINSIC :: ISO_C_BINDING
			IMPLICIT NONE
		END SUBROUTINE TwoPunctures_Run
	END INTERFACE
	
	INTERFACE
		SUBROUTINE TwoPunctures_Interpolate(x,Q) BIND(C)
			USE, INTRINSIC :: ISO_C_BINDING
			IMPLICIT NONE
			REAL, INTENT(IN)               :: x(3)
			REAL, INTENT(OUT)              :: Q(54)
		END SUBROUTINE TwoPunctures_Interpolate
	END INTERFACE
	
END MODULE TwoPunctures_mod 

#endif 
