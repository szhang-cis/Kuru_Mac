      INTERFACE
        SUBROUTINE ADDSOLT(NEQNST,NFRONT,NLASTT,NSTIFT,WORK1T)
          INCLUDE 'para_omt.f'
          INCLUDE 'addi_omt.f'
          INCLUDE 'prob_omt.f'
          INCLUDE 'inte_omt.f'
          INCLUDE 'auxl_omt.f'
          REAL*8, ALLOCATABLE, INTENT(INOUT) :: WORK1T(:)
        END SUBROUTINE
      END INTERFACE
