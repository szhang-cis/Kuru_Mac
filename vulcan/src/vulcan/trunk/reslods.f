      SUBROUTINE RESLODS(DISITS,HEADSS,IFFIXS,REFORS,TLOADS,PWORKS)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE RESIDUAL FORCES, REACTIONS & GCURN
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
      INCLUDE 'nued_om.f'
c     INCLUDE 'nuef_om.f'            ! thermal-flow  ! to be implemented
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
C
      DIMENSION DISITS(*), HEADSS(NPOINS,*), IFFIXS(*), 
     .          REFORS(*), TLOADS(NTOTVS,2)
      DIMENSION PWORKS(NPOINS,3)    ! npoins=ntotvs
C
      GCURNS=0.
      DO IDOFNS=1,NDOFNS
C$DIR NO_RECURRENCE
       DO IPOINS=1,NPOINS
        ITOTVS=(IPOINS-1)*NDOFCS+IDOFNS
c       IF(ITERME.LE.0.AND.
c    .     ITERMEF.LE.0) THEN             ! thermal or unidirec. coupled
        IF(ITERME.LE.0) THEN              ! thermal or unidirec. coupled
         IF(KINTES.EQ.0)                  ! steady-state
     .    REFORS(ITOTVS)=TLOADS(ITOTVS,1)-REFORS(ITOTVS)
         IF(KINTES.EQ.1)                  ! Euler's method
     .    REFORS(ITOTVS)=
     .    TALFAS*TLOADS(ITOTVS,1)+(1.0-TALFAS)*TLOADS(ITOTVS,2)-
     .    REFORS(ITOTVS)
         IF(KINTES.EQ.3.OR.KINTES.EQ.4)   ! Hughes' method
     .    REFORS(ITOTVS)=
     .    TLOADS(ITOTVS,1)-
     .    REFORS(ITOTVS)
        ELSE                                    ! bidirect. coupled
         IF(KINTES.EQ.0)                        ! steady-state
     .    REFORS(ITOTVS)=TLOADS(ITOTVS,1)-REFORS(ITOTVS)+
     .    PWORKS(ITOTVS,1)
         IF(KINTES.EQ.1)                        ! Euler's method
     .    REFORS(ITOTVS)=
     .    TALFAS*TLOADS(ITOTVS,1)+(1.0-TALFAS)*TLOADS(ITOTVS,2)-
     .    REFORS(ITOTVS)+
     .    TALFAS*PWORKS(ITOTVS,1)+(1.0-TALFAS)*PWORKS(ITOTVS,2)
         IF(KINTES.EQ.3.OR.KINTES.EQ.4)         ! Hughes' method
     .    REFORS(ITOTVS)=
     .    TLOADS(ITOTVS,1)-
     .    REFORS(ITOTVS)+
     .    PWORKS(ITOTVS,1)
        ENDIF        !iterme.eq.0
        GCURNS=GCURNS+REFORS(ITOTVS)*DISITS(ITOTVS)
       ENDDO
      ENDDO
C
      RETURN
      END
