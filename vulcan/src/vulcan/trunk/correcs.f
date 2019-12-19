      SUBROUTINE CORRECS(ALPHAS,DISITS,DISPRS,DISTOS,DTIMES,KINTES,
     .                   REFORS,TALFAS,TBETAS,IITERS,INDESS)
C***********************************************************************
C
C**** THIS ROUTINE UPDATES THE NODAL VARIABLES
C
C     DISIT(1:NTOTV  ) :  nodal iterative   'displacements'
C     DISPR(1:NTOTV,1) :  nodal incremental 'displacements'
C     DISPR(1:NTOTV,2) :  nodal predicted   'velocities'
C     DISPR(1:NTOTV,3) :  nodal predicted   'accelerations'
C     DISTO(1:NTOTV,1) :  nodal current     'displacements'
C     DISTO(1:NTOTV,2) :  nodal current     'velocities'
C     DISTO(1:NTOTV,3) :  nodal current     'accelerations'
C
C     INDESS=1 => correction before residual evaluation
C     INDESS=2 => correction after residual evaluation
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_oms.f'
C
      DIMENSION DISITS(*),        DISPRS(NTOTVS,*), DISTOS(NTOTVS,*),
     .          REFORS(NTOTVS,2)
C
C**** UPDATE DISPLACEMENTS, VELOCITIES AND ACCELERATIONS
C
      IF(KDYNAS.EQ.1) THEN
       CTIM1S=TALFAS/(DTIMES**2)
       CTIM2S=TBETAS/DTIMES
C
       VEDISS= 1.0D0/(TALFAS*DTIMES)
       VEVELS=-(1.0D0-TALFAS)/TALFAS
      ENDIF
C
      IF(INDESS.EQ.1) THEN       ! correction before residual evaluation
C
       DO 10 ITOTVS=1,NTOTVS
       REFORS(ITOTVS,1)  =0.0
       REFORS(ITOTVS,2)  =0.0
       DISINS         =DISPRS(ITOTVS,1)+ALPHAS*DISITS(ITOTVS)
C
       KINTEA=KINTES+1
       GO TO (1,2,3,4,5) KINTEA
C
    1  DISTOS(ITOTVS,1)=DISTOS(ITOTVS,1)+DISINS  ! steady-state
       GO TO 10
C
    2  DISTOS(ITOTVS,1)=DISTOS(ITOTVS,1)+DISINS  ! Backward Euler Method
       DISTOS(ITOTVS,2)=DISINS/DTIMES
       GO TO 10
C
    3  DISTOS(ITOTVS,1)=DISTOS(ITOTVS,1)+DISINS         ! Newmark Method
       DISTOS(ITOTVS,2)=DISPRS(ITOTVS,2)+CTIM2S*DISINS
       DISTOS(ITOTVS,3)=DISPRS(ITOTVS,3)+CTIM1S*DISINS
       GO TO 10
C
    4  DISIVS=DISPRS(ITOTVS,2)+DISITS(ITOTVS)   ! Hughes (v-form) Method
       DISTOS(ITOTVS,1)=DISTOS(ITOTVS,1)+DTIMET*DISIVS
       DISTOS(ITOTVS,2)=DISIVS
       GO TO 10
C
    5  DISTOS(ITOTVS,1)=DISTOS(ITOTVS,1)+DISINS ! Hughes (d-form) Method
       DISTOS(ITOTVS,2)=VEDISS*DISINS+VEVELS*DISPRS(ITOTVS,2)
       GO TO 10
C
   10  CONTINUE
      ENDIF
C
      IF(INDESS.EQ.2) THEN       ! correction after residual evaluation
       DO ITOTVS=1,NTOTVS
        IF(KINTES.EQ.3.AND.KPROBS.EQ.4)
     .   DISPRS(ITOTVS,2)=DISPRS(ITOTVS,2)+DISITS(ITOTVS) ! Hughes Meth.
        DISPRS(ITOTVS,1)=DISPRS(ITOTVS,1)+ALPHAS*DISITS(ITOTVS)
       ENDDO
      ENDIF
C
      RETURN
      END
