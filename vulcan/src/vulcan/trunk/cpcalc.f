      SUBROUTINE CPCALC(X1A,Y1A,Z1A,X1B,Y1B,Z1B,X1C,Y1C,Z1C,X1D,Y1D,Z1D,
     .                  X1, Y1, Z1, PX1,PY1,PZ1,
     .                  NDIME,NNODL,RCOOR,SCOOR,DISTT,DISMI,IAUXE,
     .                  RCINF,SCINF,RCSUP)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES COORDINATES OF SLAVE CONTACT POINT
C     CORRESPONDING TO THE GAUSS POINT OF THE MASTER ELEMENT
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      IAUXE=0                             ! not used !
      TOLCCT=1.0D-08                      ! better as input
      MAXITER=100
      TOLRES=1.0D-06
      TOLMAX=1.0D+04
C
      IF(NDIME.EQ.2) THEN
       DISTT=0.0D0
       DISMI=0.0D0
       RCOOR=2.0D0                        ! initial guess
       IF(NNODL.EQ.2) THEN
        ALREF=(X1A-X1B)*(X1A-X1B)+(Y1A-Y1B)*(Y1A-Y1B)
        ALREF=TOLCCT*DSQRT(ALREF)         ! reference tolerance
        DENOM=0.5D0*PX1*(Y1B-Y1A)-0.5D0*PY1*(X1B-X1A)
        IF(DABS(DENOM).LT.ALREF) RETURN   ! normal // to slave segment
        RCOOR=(0.5D0*PY1*(X1B+X1A)-PY1*X1-0.5D0*PX1*(Y1B+Y1A)+PX1*Y1)/
     .        (0.5D0*PX1*(Y1B-Y1A)-0.5D0*PY1*(X1B-X1A))
        X1AB=0.5D0*(1.0D0-RCOOR)*X1A+0.5D0*(1.0D0+RCOOR)*X1B
        Y1AB=0.5D0*(1.0D0-RCOOR)*Y1A+0.5D0*(1.0D0+RCOOR)*Y1B
        DISMI=PX1*(X1-X1AB)+PY1*(Y1-Y1AB)           ! normal gap
        DISTT=DABS(DISMI)                           ! normal gap modulus
       ENDIF                              ! nnodl.eq.2
       RETURN
      ENDIF                               ! ndime.eq.2
C
      IF(NDIME.EQ.3) THEN
       IF(NNODL.EQ.3) THEN
        DELTAR=0.0D0                      ! initialization
        DELTAS=0.0D0
        DISTT=0.0D0
        DISMI=0.0D0
        RCOOR=1.0D0                       ! initial guess
        SCOOR=1.0D0
        RCINF=0.0D0                       ! bounds
        SCINF=0.0D0
        RCSUP=RCOOR
        S=RCOOR                           ! residual evaluation
        T=SCOOR
        S1=1.0D0-S-T
        S2=S
        S3=T
        DR1=-1.0D0
        DR2= 1.0D0
        DR3= 0.0D0
        DS1=-1.0D0
        DS2= 0.0D0
        DS3= 1.0D0
        SUMX=S1*X1A+S2*X1B+S3*X1C
        SUMY=S1*Y1A+S2*Y1B+S3*Y1C
        SUMZ=S1*Z1A+S2*Z1B+S3*Z1C
C
        RESI11=PY1*(SUMX-X1)-PX1*(SUMY-Y1)    ! 3 particular cases (3pc)
        RESI21=PZ1*(SUMX-X1)-PX1*(SUMZ-Z1)
        RESI12=PX1*(SUMY-Y1)-PY1*(SUMX-X1)
        RESI22=PZ1*(SUMY-Y1)-PY1*(SUMZ-Z1)
        RESI13=PX1*(SUMZ-Z1)-PZ1*(SUMX-X1)
        RESI23=PY1*(SUMZ-Z1)-PZ1*(SUMY-Y1)
C
        DO I=1,MAXITER                    ! iterative loop (linear!)
         SUMDRX=DR1*X1A+DR2*X1B+DR3*X1C
         SUMDRY=DR1*Y1A+DR2*Y1B+DR3*Y1C
         SUMDRZ=DR1*Z1A+DR2*Z1B+DR3*Z1C
         SUMDSX=DS1*X1A+DS2*X1B+DS3*X1C
         SUMDSY=DS1*Y1A+DS2*Y1B+DS3*Y1C
         SUMDSZ=DS1*Z1A+DS2*Z1B+DS3*Z1C
C
         XJ111=-PY1*SUMDRX+PX1*SUMDRY     ! jacobian evaluation (3pc)
         XJ121=-PY1*SUMDSX+PX1*SUMDSY
         XJ211=-PZ1*SUMDRX+PX1*SUMDRZ
         XJ221=-PZ1*SUMDSX+PX1*SUMDSZ
         XJ112=-PX1*SUMDRY+PY1*SUMDRX
         XJ122=-PX1*SUMDSY+PY1*SUMDSX
         XJ212=-PZ1*SUMDRY+PY1*SUMDRZ
         XJ222=-PZ1*SUMDSY+PY1*SUMDSZ
         XJ113=-PX1*SUMDRZ+PZ1*SUMDRX
         XJ123=-PX1*SUMDSZ+PZ1*SUMDSX
         XJ213=-PY1*SUMDRZ+PZ1*SUMDRY
         XJ223=-PY1*SUMDSZ+PZ1*SUMDSY
         DETER1=XJ111*XJ221-XJ211*XJ121
         DETER2=XJ112*XJ222-XJ212*XJ122
         DETER3=XJ113*XJ223-XJ213*XJ123
C
         IMAX=1
         DETMAX=DABS(DETER1)
         IF(DABS(DETER2).GT.DETMAX) THEN
          IMAX=2
          DETMAX=DABS(DETER2)
         ENDIF
         IF(DABS(DETER3).GT.DETMAX) THEN
          IMAX=3
          DETMAX=DABS(DETER3)
         ENDIF
C
         IF(DETMAX.LT.TOLCCT) RETURN      ! normal // to slave surface
         IF(IMAX.EQ.1) THEN
          DETER=DETER1
          XJ11=XJ111
          XJ12=XJ121
          XJ21=XJ211
          XJ22=XJ221
          RESI1=RESI11
          RESI2=RESI21
         ENDIF
         IF(IMAX.EQ.2) THEN
          DETER=DETER2
          XJ11=XJ112
          XJ12=XJ122
          XJ21=XJ212
          XJ22=XJ222
          RESI1=RESI12
          RESI2=RESI22
         ENDIF
         IF(IMAX.EQ.3) THEN
          DETER=DETER3
          XJ11=XJ113
          XJ12=XJ123
          XJ21=XJ213
          XJ22=XJ223
          RESI1=RESI13
          RESI2=RESI23
         ENDIF
C
         DENOM=1.0D0/DETER
         XJI11= XJ22*DENOM                ! inverse of jacobian
         XJI22= XJ11*DENOM
         XJI21=-XJ21*DENOM
         XJI12=-XJ12*DENOM
         DELTAR=XJI11*RESI1+XJI12*RESI2   ! incremental solution
         DELTAS=XJI21*RESI1+XJI22*RESI2
         RCOOR=RCOOR+DELTAR               ! update
         SCOOR=SCOOR+DELTAS
         RCSUP=RCOOR
C
         S=RCOOR                          ! residual evaluation
         T=SCOOR
         S1=1.0D0-S-T
         S2=S
         S3=T
         DR1=-1.0D0
         DR2= 1.0D0
         DR3= 0.0D0
         DS1=-1.0D0
         DS2= 0.0D0
         DS3= 1.0D0
         SUMX=S1*X1A+S2*X1B+S3*X1C
         SUMY=S1*Y1A+S2*Y1B+S3*Y1C
         SUMZ=S1*Z1A+S2*Z1B+S3*Z1C
C
         RESI11=PY1*(SUMX-X1)-PX1*(SUMY-Y1)
         RESI21=PZ1*(SUMX-X1)-PX1*(SUMZ-Z1)
         RESI12=PX1*(SUMY-Y1)-PY1*(SUMX-X1)
         RESI22=PZ1*(SUMY-Y1)-PY1*(SUMZ-Z1)
         RESI13=PX1*(SUMZ-Z1)-PZ1*(SUMX-X1)
         RESI23=PY1*(SUMZ-Z1)-PZ1*(SUMY-Y1)
C
         IF(IMAX.EQ.1) THEN
          RESI1=RESI11
          RESI2=RESI21
         ENDIF
         IF(IMAX.EQ.2) THEN
          RESI1=RESI12
          RESI2=RESI22
         ENDIF
         IF(IMAX.EQ.3) THEN
          RESI1=RESI13
          RESI2=RESI23
         ENDIF
         RESIT=RESI1*RESI1+RESI2*RESI2    ! convergence criterion
         IF(DABS(RESIT).LT.TOLRES) GO TO 100
c        IF(I.EQ.MAXITER) IAUXE=1
         IF(I.EQ.MAXITER) THEN            ! intersection not found: skip
          RCOOR=1.0D0
          SCOOR=1.0D0
          RETURN
         ENDIF
        ENDDO
  100   CONTINUE
        X1ABC=(1.0D0-RCOOR-SCOOR)*X1A+RCOOR*X1B+SCOOR*X1C
        Y1ABC=(1.0D0-RCOOR-SCOOR)*Y1A+RCOOR*Y1B+SCOOR*Y1C
        Z1ABC=(1.0D0-RCOOR-SCOOR)*Z1A+RCOOR*Z1B+SCOOR*Z1C
        DISMI=PX1*(X1-X1ABC)+PY1*(Y1-Y1ABC)+        ! normal gap
     .        PZ1*(Z1-Z1ABC)
        DISTT=DABS(DISMI)                           ! normal gap modulus
       ENDIF                              ! nnodl.eq.3
C
       IF(NNODL.EQ.4) THEN                ! Newton-Raphson algorithm
        DELTAR=0.0D0                      ! initialization
        DELTAS=0.0D0
        DISTT=0.0D0
        DISMI=0.0D0
        RCOOR=0.0D0                       ! initial guess
        SCOOR=0.0D0
        RCINF=-1.0D0                      ! bounds
        SCINF=-1.0D0
        RCSUP=0.0D0
        S=RCOOR                           ! residual evaluation
        T=SCOOR
        ST=S*T
        S1=(1.0D0-T-S+ST)*0.25D0
        S2=(1.0D0-T+S-ST)*0.25D0
        S3=(1.0D0+T+S+ST)*0.25D0
        S4=(1.0D0+T-S-ST)*0.25D0
        DR1=(-1.0D0+T)*0.25D0
        DR2=(+1.0D0-T)*0.25D0
        DR3=(+1.0D0+T)*0.25D0
        DR4=(-1.0D0-T)*0.25D0
        DS1=(-1.0D0+S)*0.25D0
        DS2=(-1.0D0-S)*0.25D0
        DS3=(+1.0D0+S)*0.25D0
        DS4=(+1.0D0-S)*0.25D0
        SUMX=S1*X1A+S2*X1B+S3*X1C+S4*X1D
        SUMY=S1*Y1A+S2*Y1B+S3*Y1C+S4*Y1D
        SUMZ=S1*Z1A+S2*Z1B+S3*Z1C+S4*Z1D
C
        RESI11=PY1*(SUMX-X1)-PX1*(SUMY-Y1)    ! 3 particular cases (3pc)
        RESI21=PZ1*(SUMX-X1)-PX1*(SUMZ-Z1)
        RESI12=PX1*(SUMY-Y1)-PY1*(SUMX-X1)
        RESI22=PZ1*(SUMY-Y1)-PY1*(SUMZ-Z1)
        RESI13=PX1*(SUMZ-Z1)-PZ1*(SUMX-X1)
        RESI23=PY1*(SUMZ-Z1)-PZ1*(SUMY-Y1)
C
        DO I=1,MAXITER                    ! iterative loop
         SUMDRX=DR1*X1A+DR2*X1B+DR3*X1C+DR4*X1D
         SUMDRY=DR1*Y1A+DR2*Y1B+DR3*Y1C+DR4*Y1D
         SUMDRZ=DR1*Z1A+DR2*Z1B+DR3*Z1C+DR4*Z1D
         SUMDSX=DS1*X1A+DS2*X1B+DS3*X1C+DS4*X1D
         SUMDSY=DS1*Y1A+DS2*Y1B+DS3*Y1C+DS4*Y1D
         SUMDSZ=DS1*Z1A+DS2*Z1B+DS3*Z1C+DS4*Z1D
C
         XJ111=-PY1*SUMDRX+PX1*SUMDRY     ! jacobian evaluation (3pc)
         XJ121=-PY1*SUMDSX+PX1*SUMDSY
         XJ211=-PZ1*SUMDRX+PX1*SUMDRZ
         XJ221=-PZ1*SUMDSX+PX1*SUMDSZ
         XJ112=-PX1*SUMDRY+PY1*SUMDRX
         XJ122=-PX1*SUMDSY+PY1*SUMDSX
         XJ212=-PZ1*SUMDRY+PY1*SUMDRZ
         XJ222=-PZ1*SUMDSY+PY1*SUMDSZ
         XJ113=-PX1*SUMDRZ+PZ1*SUMDRX
         XJ123=-PX1*SUMDSZ+PZ1*SUMDSX
         XJ213=-PY1*SUMDRZ+PZ1*SUMDRY
         XJ223=-PY1*SUMDSZ+PZ1*SUMDSY
         DETER1=XJ111*XJ221-XJ211*XJ121
         DETER2=XJ112*XJ222-XJ212*XJ122
         DETER3=XJ113*XJ223-XJ213*XJ123
C
         IMAX=1
         DETMAX=DABS(DETER1)
         IF(DABS(DETER2).GT.DETMAX) THEN
          IMAX=2
          DETMAX=DABS(DETER2)
         ENDIF
         IF(DABS(DETER3).GT.DETMAX) THEN
          IMAX=3
          DETMAX=DABS(DETER3)
         ENDIF
C
         IF(DETMAX.LT.TOLCCT) THEN        ! normal // to slave surface
          RCOOR=2.0D0
          SCOOR=2.0D0
          RETURN
         ENDIF
         IF(IMAX.EQ.1) THEN
          DETER=DETER1
          XJ11=XJ111
          XJ12=XJ121
          XJ21=XJ211
          XJ22=XJ221
          RESI1=RESI11
          RESI2=RESI21
         ENDIF
         IF(IMAX.EQ.2) THEN
          DETER=DETER2
          XJ11=XJ112
          XJ12=XJ122
          XJ21=XJ212
          XJ22=XJ222
          RESI1=RESI12
          RESI2=RESI22
         ENDIF
         IF(IMAX.EQ.3) THEN
          DETER=DETER3
          XJ11=XJ113
          XJ12=XJ123
          XJ21=XJ213
          XJ22=XJ223
          RESI1=RESI13
          RESI2=RESI23
         ENDIF
C
         DENOM=1.0D0/DETER
         XJI11= XJ22*DENOM                ! inverse of jacobian
         XJI22= XJ11*DENOM
         XJI21=-XJ21*DENOM
         XJI12=-XJ12*DENOM
         DELTAR=XJI11*RESI1+XJI12*RESI2   ! incremental solution
         DELTAS=XJI21*RESI1+XJI22*RESI2
         RCOOR=RCOOR+DELTAR               ! update
         SCOOR=SCOOR+DELTAS
C
         IF(DABS(RCOOR).GT.TOLMAX.OR.
     .      DABS(SCOOR).GT.TOLMAX) THEN   ! intersection not found: skip
          RCOOR=2.0D0
          SCOOR=2.0D0
          RETURN
         ENDIF
C
         S=RCOOR                          ! residual evaluation
         T=SCOOR
         ST=S*T
         S1=(1.0D0-T-S+ST)*0.25D0
         S2=(1.0D0-T+S-ST)*0.25D0
         S3=(1.0D0+T+S+ST)*0.25D0
         S4=(1.0D0+T-S-ST)*0.25D0
         DR1=(-1.0D0+T)*0.25D0
         DR2=(+1.0D0-T)*0.25D0
         DR3=(+1.0D0+T)*0.25D0
         DR4=(-1.0D0-T)*0.25D0
         DS1=(-1.0D0+S)*0.25D0
         DS2=(-1.0D0-S)*0.25D0
         DS3=(+1.0D0+S)*0.25D0
         DS4=(+1.0D0-S)*0.25D0
         SUMX=S1*X1A+S2*X1B+S3*X1C+S4*X1D
         SUMY=S1*Y1A+S2*Y1B+S3*Y1C+S4*Y1D
         SUMZ=S1*Z1A+S2*Z1B+S3*Z1C+S4*Z1D
C
         RESI11=PY1*(SUMX-X1)-PX1*(SUMY-Y1)
         RESI21=PZ1*(SUMX-X1)-PX1*(SUMZ-Z1)
         RESI12=PX1*(SUMY-Y1)-PY1*(SUMX-X1)
         RESI22=PZ1*(SUMY-Y1)-PY1*(SUMZ-Z1)
         RESI13=PX1*(SUMZ-Z1)-PZ1*(SUMX-X1)
         RESI23=PY1*(SUMZ-Z1)-PZ1*(SUMY-Y1)
C
         IF(IMAX.EQ.1) THEN
          RESI1=RESI11
          RESI2=RESI21
         ENDIF
         IF(IMAX.EQ.2) THEN
          RESI1=RESI12
          RESI2=RESI22
         ENDIF
         IF(IMAX.EQ.3) THEN
          RESI1=RESI13
          RESI2=RESI23
         ENDIF
         RESIT=RESI1*RESI1+RESI2*RESI2    ! convergence criterion
         IF(DABS(RESIT).LT.TOLRES) GO TO 101
c        IF(I.EQ.MAXITER) IAUXE=1
         IF(I.EQ.MAXITER) THEN            ! intersection not found: skip
          RCOOR=2.0D0
          SCOOR=2.0D0
          RETURN
         ENDIF
        ENDDO
  101   CONTINUE
        X1ABCD=0.25D0*(1.0D0-RCOOR)*(1.0D0-SCOOR)*X1A+
     .         0.25D0*(1.0D0+RCOOR)*(1.0D0-SCOOR)*X1B+
     .         0.25D0*(1.0D0+RCOOR)*(1.0D0+SCOOR)*X1C+
     .         0.25D0*(1.0D0-RCOOR)*(1.0D0+SCOOR)*X1D
        Y1ABCD=0.25D0*(1.0D0-RCOOR)*(1.0D0-SCOOR)*Y1A+
     .         0.25D0*(1.0D0+RCOOR)*(1.0D0-SCOOR)*Y1B+
     .         0.25D0*(1.0D0+RCOOR)*(1.0D0+SCOOR)*Y1C+
     .         0.25D0*(1.0D0-RCOOR)*(1.0D0+SCOOR)*Y1D
        Z1ABCD=0.25D0*(1.0D0-RCOOR)*(1.0D0-SCOOR)*Z1A+
     .         0.25D0*(1.0D0+RCOOR)*(1.0D0-SCOOR)*Z1B+
     .         0.25D0*(1.0D0+RCOOR)*(1.0D0+SCOOR)*Z1C+
     .         0.25D0*(1.0D0-RCOOR)*(1.0D0+SCOOR)*Z1D
        DISMI=PX1*(X1-X1ABCD)+PY1*(Y1-Y1ABCD)+      ! normal gap
     .        PZ1*(Z1-Z1ABCD)
        DISTT=DABS(DISMI)                           ! normal gap modulus
       ENDIF                              ! nnodl.eq.4
       RETURN
      ENDIF                               ! ndime.eq.3
C
      END
