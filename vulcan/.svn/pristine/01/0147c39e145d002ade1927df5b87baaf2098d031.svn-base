      SUBROUTINE PREDICS(DISITS,DISPRS,DISTOS)
C***********************************************************************
C
C**** THIS ROUTINE MAKES A PREDICTION ON NODAL VARIABLES
C
C          DISITT(1:NTOTVT  ) :  nodal iterative   'displacements'
C          DISPRT(1:NTOTVT,1) :  nodal incremental 'displacements'
C          DISPRT(1:NTOTVT,2) :  nodal predicted   'velocities'
C          DISPRT(1:NTOTVT,3) :  nodal predicted   'accelerations'
C          DISTOT(1:NTOTVT,1) :  nodal current     'displacements'
C          DISTOT(1:NTOTVT,2) :  nodal current     'velocities'
C          DISTOT(1:NTOTVT,3) :  nodal current     'accelerations'
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
C
      DIMENSION DISITS(*), DISPRS(NTOTVS,*), DISTOS(NTOTVS,*)
C
C**** INITIALISE ITERATIVE DISPLACEMENTS
C
      SCALES=0.0D+00
cc    IF(FACPRS.NE.0.0) SCALES=FACTOS/FACPRS
C
      IF(LACCES.NE.0.AND.ISTEPS.GT.1.AND.SCALES.NE.0.0) THEN      
        DO 10 ITOTVS=1,NTOTVS
   10   DISITS(ITOTVS)=DISPRS(ITOTVS,1)*SCALES
      ELSE
        DO 20 ITOTVS=1,NTOTVS
   20   DISITS(ITOTVS)=0.0
      ENDIF
C
C**** INITIALISE INCREMENTAL DISPLACEMENTS  
C
      DO 30 ITOTVS=1,NTOTVS
      DISPRS(ITOTVS,1)=0.0
C
      IF(KINTES.EQ.1) THEN   ! PREDICT VELOC. (EULER'S METHOD)
       DISPRS(ITOTVS,2)= DISTOS(ITOTVS,2)
      ENDIF
C
      IF(KINTES.EQ.2) THEN   ! PREDICT VELOC. & ACCEL.(NEWMARK'S METHOD)
       DISPRS(ITOTVS,2)=  (1.-   TBETAS)*DISTOS(ITOTVS,2)
     .                   +(1.-.5*TBETAS)*DISTOS(ITOTVS,3)*DTIMES
       DISPRS(ITOTVS,3)=        -TALFAS *DISTOS(ITOTVS,2)/DTIMES
     .                   +(1.-.5*TALFAS)*DISTOS(ITOTVS,3)
      ENDIF
C
      IF(KINTES.EQ.3)
     . CALL RUNENDS('HIGHES-V_FORM NOT IMPLEMENTED')
C
      IF(KINTES.EQ.4) THEN   ! PREDICT VELOC. (HUGHES-D_FORM METHOD)
       DISPRS(ITOTVS,2)=DISTOS(ITOTVS,2)
      ENDIF
C
   30 CONTINUE    
C
      RETURN
      END
