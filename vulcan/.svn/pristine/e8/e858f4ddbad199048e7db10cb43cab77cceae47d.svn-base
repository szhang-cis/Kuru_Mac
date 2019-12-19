      SUBROUTINE PREDICT(DISITT,DISPRT,DISTOT)
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
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION DISITT(*), DISPRT(NTOTVT,*), DISTOT(NTOTVT,*)
C
C**** INITIALISE ITERATIVE DISPLACEMENTS
C
      SCALET=0.0D+00
cc      IF(FACPRT.NE.0.0) SCALET=FACTOT/FACPRT
C
      IF(LACCET.NE.0.AND.ISTEPT.GT.1.AND.SCALET.NE.0.0) THEN      
        DO 10 ITOTVT=1,NTOTVT
   10   DISITT(ITOTVT)=DISPRT(ITOTVT,1)*SCALET
      ELSE
        DO 20 ITOTVT=1,NTOTVT
   20   DISITT(ITOTVT)=0.0
      ENDIF
C
C**** INITIALISE INCREMENTAL DISPLACEMENTS  
C
      DO 30 ITOTVT=1,NTOTVT
      DISPRT(ITOTVT,1)=0.0
C
      IF(KINTET.EQ.1) THEN   ! PREDICT VELOC. (EULER'S METHOD)
       DISPRT(ITOTVT,2)= DISTOT(ITOTVT,2)
      ENDIF
C
      IF(KINTET.EQ.2) THEN   ! PREDICT VELOC. & ACCEL.(NEWMARK'S METHOD)
       DISPRT(ITOTVT,2)=  (1.-   TBETAT)*DISTOT(ITOTVT,2)
     .                   +(1.-.5*TBETAT)*DISTOT(ITOTVT,3)*DTIMET
       DISPRT(ITOTVT,3)=        -TALFAT *DISTOT(ITOTVT,2)/DTIMET
     .                   +(1.-.5*TALFAT)*DISTOT(ITOTVT,3)
      ENDIF
C
      IF(KINTET.EQ.3)
     . CALL RUNENDT('HIGHES-V_FORM NOT IMPLEMENTED')
C
      IF(KINTET.EQ.4) THEN   ! PREDICT VELOC. (HUGHES-D_FORM METHOD)
       DISPRT(ITOTVT,2)=DISTOT(ITOTVT,2)
      ENDIF
C
   30 CONTINUE    
C
      RETURN
      END
