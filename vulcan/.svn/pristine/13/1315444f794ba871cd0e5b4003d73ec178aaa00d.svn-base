      SUBROUTINE ASELM3S(MATNOS,PROELS,PROPSS,CSTIFS,ESTIFS,WSTIFS)
C***********************************************************************
C
C**** THIS ROUTINE ASSEMBLES ELEMENTAL "EFFECTIVE" MATRIX  FOR
C     UNCOUPLED PROBLEMS
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
C
      DIMENSION MATNOS(*),             PROELS(NPRELS,*),
     .          PROPSS(NPROPS,*)
      DIMENSION CSTIFS(NEVABS,NEVABS), ESTIFS(NKOVAS),
     .          WSTIFS(NKOVAS)
C
C**** ASSEMBLE ELEMENTAL MATRIX
C
      IF(KDYNAS.EQ.1) THEN       ! Euler Method
       CTIM2S=1.0D00/DTIMES
ctm    CTIM3S=1.0D00
       CTIM3S=TALFAS
       IF(KINTES.EQ.3) THEN      ! Hughes (v-form) Method
        CTIM2S=1.0D00
        CTIM3S=TALFAS*DTIMES
       ENDIF
       IF(KINTES.EQ.4) THEN      ! Hughes (d-form) Method
        CTIM2S=1.0D00/(TALFAS*DTIMES)
        CTIM3S=1.0D00
       ENDIF
      ENDIF
C
C**** LOOP OVER ELEMENTS
C
      DO 1000 IELEMS=1,NELEMS
      LGRUPS=MATNOS(IELEMS)
      LMATSS=INT(PROELS(1,LGRUPS))
C
C**** READ MATRICES FROM DATA BASE
C 
                      CALL DATBASS(ESTIFS,    7,    2)
      IF(KDYNAS.EQ.1) CALL DATBASS(WSTIFS,    8,    2)
C
C**** CHECK THE CORRECTNESS OF THE CONSTITUENT ELEMENTAL MATRICES
C
      ICHEKS=0
      IF(ICHEKS.EQ.1) THEN
        IF(KPROBS.GE.1)
     .    CALL TESTM1(ESTIFS,IELEMS,NEVABS)
        IF(KDYNAS.NE.0)
     .    CALL TESTM1(WSTIFS,IELEMS,NEVABS)
      ENDIF
C
C**** INITIALISE ASSEMBLED ELEMENTAL MATRIX
C
      DO IEVABS=1,NEVABS
       DO JEVABS=1,NEVABS
        CSTIFS(IEVABS,JEVABS)=0.0
       ENDDO
      ENDDO
C
      IF(KSYMMS.EQ.0) THEN          ! unsymmetric case
C
       IF(KDYNAS.EQ.0) THEN         ! steady-state case
        IKOVAS=0
        DO IEVABS=1,NEVABS
         DO JEVABS=1,NEVABS
          IKOVAS=IKOVAS+1
          CSTIFS(IEVABS,JEVABS)=ESTIFS(IKOVAS)
         ENDDO
        ENDDO
       ELSE                         ! transient case
        IKOVAS=0
        DO IEVABS=1,NEVABS
         DO JEVABS=1,NEVABS
          IKOVAS=IKOVAS+1 
          CSTIFS(IEVABS,JEVABS)=CTIM2S*WSTIFS(IKOVAS)+
     .                          CTIM3S*ESTIFS(IKOVAS)
         ENDDO
        ENDDO
       ENDIF
C
      ELSE                          ! symmetric case
C
       IF(KDYNAS.EQ.0) THEN         ! steady-state case
        IKOVAS=0
        DO IEVABS=1,NEVABS
         DO JEVABS=IEVABS,NEVABS
          IKOVAS=IKOVAS+1
          CSTIFS(IEVABS,JEVABS)=ESTIFS(IKOVAS)
          CSTIFS(JEVABS,IEVABS)=CSTIFS(IEVABS,JEVABS)
         ENDDO
        ENDDO
       ELSE                         ! transient case
        IKOVAS=0
        DO IEVABS=1,NEVABS
         DO JEVABS=IEVABS,NEVABS
          IKOVAS=IKOVAS+1 
          CSTIFS(IEVABS,JEVABS)=CTIM2S*WSTIFS(IKOVAS)+
     .                          CTIM3S*ESTIFS(IKOVAS)
          CSTIFS(JEVABS,IEVABS)=CSTIFS(IEVABS,JEVABS)
         ENDDO
        ENDDO
       ENDIF
C
      ENDIF
C
C**** WRITE CSTIF TO DATA BASE
C 
      CALL DATBASS(CSTIFS,    6,    1)
C
 1000 CONTINUE
C
      RETURN
      END
