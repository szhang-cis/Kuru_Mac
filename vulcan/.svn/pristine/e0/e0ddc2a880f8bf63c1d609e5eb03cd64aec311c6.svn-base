      SUBROUTINE ASELM3T(MATNOT,PROELT,PROPST,CSTIFT,ESTIFT,WSTIFT)
C***********************************************************************
C
C**** THIS ROUTINE ASSEMBLES ELEMENTAL "EFFECTIVE" MATRIX  FOR
C     UNCOUPLED PROBLEMS
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION MATNOT(*),             PROELT(NPRELT,*),
     .          PROPST(NPROPT,*)
      DIMENSION CSTIFT(NEVABT,NEVABT), ESTIFT(NKOVAT),
     .          WSTIFT(NKOVAT)
C
C**** ASSEMBLE ELEMENTAL MATRIX
C
      IF(KDYNAT.EQ.1) THEN       ! Euler Method
       CTIM2T=1.0D00/DTIMET
ctm    CTIM3T=1.0D00
       CTIM3T= TALFAT
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
C**** LOOP OVER ELEMENTS
C
      DO 1000 IELEMT=1,NELEMT
      LGRUPT=MATNOT(IELEMT)
      LMATST=INT(PROELT(1,LGRUPT))
C
C**** READ MATRICES FROM DATA BASE
C 
                      CALL DATBAST(ESTIFT,    7,    2)
      IF(KDYNAT.EQ.1) CALL DATBAST(WSTIFT,    8,    2)
C
C**** CHECK THE CORRECTNESS OF THE CONSTITUENT ELEMENTAL MATRICES
C
      ICHEKT=0
      IF(ICHEKT.EQ.1) THEN
        IF(KPROBT.GE.1)
     .    CALL TESTM1(ESTIFT,IELEMT,NEVABT)
        IF(KDYNAT.NE.0)
     .    CALL TESTM1(WSTIFT,IELEMT,NEVABT)
      ENDIF
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
C**** WRITE CSTIF TO DATA BASE
C 
      CALL DATBAST(CSTIFT,    6,    1)
C
 1000 CONTINUE
C
      RETURN
      END
