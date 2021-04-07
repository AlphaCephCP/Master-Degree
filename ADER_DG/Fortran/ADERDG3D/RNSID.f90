#ifdef RNSID 
!
! This is a Fortran adapter to the RNSID C++ code.
! It makes use of the plain C RNSIDadapter.c.
!
! This is only needed when you want to use RNSID from your
! Fortran code.
!
    
MODULE RNSID_mod  
        ! 
        IMPLICIT NONE
        PUBLIC 
        ! 
        INTERFACE
                SUBROUTINE RNSID_Setup() BIND(C)
                        USE, INTRINSIC :: ISO_C_BINDING
                        IMPLICIT NONE
                END SUBROUTINE RNSID_Setup
        END INTERFACE

        INTERFACE
                SUBROUTINE RNSID_AmendParameters() BIND(C)
                        USE, INTRINSIC :: ISO_C_BINDING
                        IMPLICIT NONE
                END SUBROUTINE RNSID_AmendParameters
        END INTERFACE

        INTERFACE
                SUBROUTINE RNSID_Run() BIND(C)
                        USE, INTRINSIC :: ISO_C_BINDING
                        IMPLICIT NONE
                END SUBROUTINE RNSID_Run
        END INTERFACE

        INTERFACE
                SUBROUTINE RNSID_Interpolate(x,Q) BIND(C)
                        USE, INTRINSIC :: ISO_C_BINDING
                        IMPLICIT NONE
                        REAL, INTENT(IN)               :: x(3)
                        REAL, INTENT(OUT)              :: Q(70)   ! Tune this to your length of the conserved quantities
                END SUBROUTINE RNSID_Interpolate
        END INTERFACE
        ! 
END MODULE RNSID_mod  
    

#endif 

