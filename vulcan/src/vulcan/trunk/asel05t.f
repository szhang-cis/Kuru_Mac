      SUBROUTINE ASEL05T(CSTIFT,ESTIFT,WSTIFT,CSTI1T)
C***********************************************************************
C
C**** THIS ROUTINE ASSEMBLES ELEMENTAL "EFFECTIVE" MATRIX  FOR
C     UNCOUPLED PROBLEMS
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION CSTIFT(NEVABT,NEVABT), ESTIFT(NKOVAT),
     .          WSTIFT(NKOVAT)
      DIMENSION CSTI1T(NEVABT,NEVABT)
C
      IF(NMEMO6.EQ.0) RETURN
C
C**** ASSEMBLE ELEMENTAL MATRIX
C
      IF(KDYNAT.EQ.1) THEN       ! Euler Method
       CTIM2T=1.0D00/DTIMET
ctm    CTIM3T=1.0D00
       CTIM3T=TALFAT
       IF(KINTET.EQ.3) THEN      ! Hughes (v-form) Method
        CTIM2T=1.0D00
        CTIM3T=TALFAT*DTIMET
       ENDIF
       IF(KINTET.EQ.4) THEN      ! Hughes (d-form) Method
        CTIM2T=1.0D00/(TALFAT*DTIMET)
        CTIM3T=1.0D00
       ENDIF
      ENDIF
C
C**** CHECK THE CORRECTNESS OF THE CONSTITUENT ELEMENTAL MATRICES
C
      ICHEKT=0
      IF(ICHEKT.EQ.1) THEN
       IF(KPROBT.GE.1)
     .  CALL TESTM1(ESTIFT,IELEMT,NEVABT)
       IF(KDYNAT.NE.0)
     .  CALL TESTM1(WSTIFT,IELEMT,NEVABT)
      ENDIF
C
      IF(NMEMO7.EQ.0) THEN
C
C**** INITIALISE ASSEMBLED ELEMENTAL MATRIX
C
      DO IEVABT=1,NEVABT
       DO JEVABT=1,NEVABT
        CSTIFT(IEVABT,JEVABT)=0.0
       ENDDO
      ENDDO
C
      IF(KSYMMT.EQ.0) THEN          ! unsymmetric case
C
       IF(KDYNAT.EQ.0) THEN         ! steady-state case
        IKOVAT=0
        DO IEVABT=1,NEVABT
         DO JEVABT=1,NEVABT
          IKOVAT=IKOVAT+1
          CSTIFT(IEVABT,JEVABT)=ESTIFT(IKOVAT)
         ENDDO
        ENDDO
       ELSE                         ! transient case
        IKOVAT=0
        DO IEVABT=1,NEVABT
         DO JEVABT=1,NEVABT
          IKOVAT=IKOVAT+1 
          CSTIFT(IEVABT,JEVABT)=CTIM2T*WSTIFT(IKOVAT)+
     .                          CTIM3T*ESTIFT(IKOVAT)
         ENDDO
        ENDDO
       ENDIF
C
      ELSE                          ! symmetric case
C
       IF(KDYNAT.EQ.0) THEN         ! steady-state case
        IKOVAT=0
        DO IEVABT=1,NEVABT
         DO JEVABT=IEVABT,NEVABT
          IKOVAT=IKOVAT+1
          CSTIFT(IEVABT,JEVABT)=ESTIFT(IKOVAT)
          CSTIFT(JEVABT,IEVABT)=CSTIFT(IEVABT,JEVABT)
         ENDDO
        ENDDO
       ELSE                         ! transient case
        IKOVAT=0
        DO IEVABT=1,NEVABT
         DO JEVABT=IEVABT,NEVABT
          IKOVAT=IKOVAT+1 
          CSTIFT(IEVABT,JEVABT)=CTIM2T*WSTIFT(IKOVAT)+
     .                          CTIM3T*ESTIFT(IKOVAT)
          CSTIFT(JEVABT,IEVABT)=CSTIFT(IEVABT,JEVABT)
         ENDDO
        ENDDO
       ENDIF
C
      ENDIF
C
      ELSE                  ! nmemo7=1
C
C**** INITIALISE ASSEMBLED ELEMENTAL MATRIX
C
      DO IEVABT=1,NEVABT
       DO JEVABT=1,NEVABT
        CSTI1T(IEVABT,JEVABT)=0.0
       ENDDO
      ENDDO
C
      IF(KSYMMT.EQ.0) THEN          ! unsymmetric case
C
       IF(KDYNAT.EQ.0) THEN         ! steady-state case
        IKOVAT=0
        DO IEVABT=1,NEVABT
         DO JEVABT=1,NEVABT
          IKOVAT=IKOVAT+1
          CSTI1T(IEVABT,JEVABT)=ESTIFT(IKOVAT)
         ENDDO
        ENDDO
       ELSE                         ! transient case
        IKOVAT=0
        DO IEVABT=1,NEVABT
         DO JEVABT=1,NEVABT
          IKOVAT=IKOVAT+1 
          CSTI1T(IEVABT,JEVABT)=CTIM2T*WSTIFT(IKOVAT)+
     .                          CTIM3T*ESTIFT(IKOVAT)
         ENDDO
        ENDDO
       ENDIF
C
      ELSE                          ! symmetric case
C
       IF(KDYNAT.EQ.0) THEN         ! steady-state case
        IKOVAT=0
        DO IEVABT=1,NEVABT
         DO JEVABT=IEVABT,NEVABT
          IKOVAT=IKOVAT+1
          CSTI1T(IEVABT,JEVABT)=ESTIFT(IKOVAT)
          CSTI1T(JEVABT,IEVABT)=CSTI1T(IEVABT,JEVABT)
         ENDDO
        ENDDO
       ELSE                         ! transient case
        IKOVAT=0
        DO IEVABT=1,NEVABT
         DO JEVABT=IEVABT,NEVABT
          IKOVAT=IKOVAT+1 
          CSTI1T(IEVABT,JEVABT)=CTIM2T*WSTIFT(IKOVAT)+
     .                          CTIM3T*ESTIFT(IKOVAT)
          CSTI1T(JEVABT,IEVABT)=CSTI1T(IEVABT,JEVABT)
         ENDDO
        ENDDO
       ENDIF
C
      ENDIF
C
      ENDIF             ! nmemo7.eq.0
C
      RETURN
      END
