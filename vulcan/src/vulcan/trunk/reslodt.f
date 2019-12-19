      SUBROUTINE RESLODT(DISITT,HEADST,IFFIXT,REFORT,TLOADT,PWORKT)
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
c     INCLUDE 'nuef_om.f'            ! thermal-flow  ! to be implemented
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION DISITT(*), HEADST(NPOINT,*), IFFIXT(*), 
     .          REFORT(*), TLOADT(NTOTVT,2)
      DIMENSION PWORKT(NPOINT,3)    ! npoint=ntotvt
C
      GCURNT=0.
      DO IDOFNT=1,NDOFNT
C$DIR NO_RECURRENCE
       DO IPOINT=1,NPOINT
        ITOTVT=(IPOINT-1)*NDOFCT+IDOFNT
c       IF(ITERME.LE.0.AND.
c    .     ITERMEF.LE.0) THEN             ! thermal or unidirec. coupled
        IF(ITERME.LE.0) THEN              ! thermal or unidirec. coupled
         IF(KINTET.EQ.0)                  ! steady-state
     .    REFORT(ITOTVT)=TLOADT(ITOTVT,1)-REFORT(ITOTVT)
         IF(KINTET.EQ.1)                  ! Euler's method
     .    REFORT(ITOTVT)=
     .    TALFAT*TLOADT(ITOTVT,1)+(1.0-TALFAT)*TLOADT(ITOTVT,2)-
     .    REFORT(ITOTVT)
         IF(KINTET.EQ.3.OR.KINTET.EQ.4)   ! Hughes' method
     .    REFORT(ITOTVT)=
     .    TLOADT(ITOTVT,1)-
     .    REFORT(ITOTVT)
        ELSE                                    ! bidirect. coupled
         IF(KINTET.EQ.0)                        ! steady-state
     .    REFORT(ITOTVT)=TLOADT(ITOTVT,1)-REFORT(ITOTVT)+
     .    PWORKT(ITOTVT,1)
         IF(KINTET.EQ.1)                        ! Euler's method
     .    REFORT(ITOTVT)=
     .    TALFAT*TLOADT(ITOTVT,1)+(1.0-TALFAT)*TLOADT(ITOTVT,2)-
     .    REFORT(ITOTVT)+
     .    TALFAT*PWORKT(ITOTVT,1)+(1.0-TALFAT)*PWORKT(ITOTVT,2)
         IF(KINTET.EQ.3.OR.KINTET.EQ.4)         ! Hughes' method
     .    REFORT(ITOTVT)=
     .    TLOADT(ITOTVT,1)-
     .    REFORT(ITOTVT)+
     .    PWORKT(ITOTVT,1)
        ENDIF        !iterme.eq.0
        GCURNT=GCURNT+REFORT(ITOTVT)*DISITT(ITOTVT)
       ENDDO
      ENDDO
C
      RETURN
      END
