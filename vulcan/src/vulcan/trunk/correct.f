      SUBROUTINE CORRECT(ALPHAT,DISITT,DISPRT,DISTOT,DTIMET,KINTET,
     .                   REFORT,TALFAT,TBETAT,IITERT,INDEST)
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
C     INDEST=1 => correction before residual evaluation
C     INDEST=2 => correction after residual evaluation
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_omt.f'
C
      DIMENSION DISITT(*),        DISPRT(NTOTVT,*), DISTOT(NTOTVT,*),
     .          REFORT(NTOTVT,2)
C
      IF(KDYNAT.EQ.1) THEN
       CTIM1T=TALFAT/(DTIMET**2)
       CTIM2T=TBETAT/DTIMET
C
       VEDIST= 1.0D0/(TALFAT*DTIMET)
       VEVELT=-(1.0D0-TALFAT)/TALFAT
      ENDIF
C
      IF(INDEST.EQ.1) THEN       ! correction before residual evaluation
C
C**** UPDATE DISPLACEMENTS, VELOCITIES AND ACCELERATIONS
C
       DO 10 ITOTVT=1,NTOTVT
       REFORT(ITOTVT,1)  =0.0
       REFORT(ITOTVT,2)  =0.0
       DISINT         =DISPRT(ITOTVT,1)+ALPHAT*DISITT(ITOTVT)
C
       KINTEA=KINTET+1
       GO TO (1,2,3,4,5) KINTEA
C
    1  DISTOT(ITOTVT,1)=DISTOT(ITOTVT,1)+DISINT  ! steady-state
       GO TO 10
C
    2  DISTOT(ITOTVT,1)=DISTOT(ITOTVT,1)+DISINT  ! Backward Euler Method
       DISTOT(ITOTVT,2)=DISINT/DTIMET
       GO TO 10
C
    3  DISTOT(ITOTVT,1)=DISTOT(ITOTVT,1)+DISINT         ! Newmark Method
       DISTOT(ITOTVT,2)=DISPRT(ITOTVT,2)+CTIM2T*DISINT
       DISTOT(ITOTVT,3)=DISPRT(ITOTVT,3)+CTIM1T*DISINT
       GO TO 10
C
    4  DISIVT=DISPRT(ITOTVT,2)+DISITT(ITOTVT)   ! Hughes (v-form) Method
       DISTOT(ITOTVT,1)=DISTOT(ITOTVT,1)+DTIMET*DISIVT
       DISTOT(ITOTVT,2)=DISIVT
       GO TO 10
C
    5  DISTOT(ITOTVT,1)=DISTOT(ITOTVT,1)+DISINT ! Hughes (d-form) Method
       DISTOT(ITOTVT,2)=VEDIST*DISINT+VEVELT*DISPRT(ITOTVT,2)
       GO TO 10
C
   10  CONTINUE
      ENDIF
C
      IF(INDEST.EQ.2) THEN       ! correction after residual evaluation
       DO ITOTVT=1,NTOTVT
        IF(KINTET.EQ.3.AND.KPROBT.EQ.4)
     .   DISPRT(ITOTVT,2)=DISPRT(ITOTVT,2)+DISITT(ITOTVT) ! Hughes Meth.
        DISPRT(ITOTVT,1)=DISPRT(ITOTVT,1)+ALPHAT*DISITT(ITOTVT)
       ENDDO
      ENDIF
C
      RETURN
      END
