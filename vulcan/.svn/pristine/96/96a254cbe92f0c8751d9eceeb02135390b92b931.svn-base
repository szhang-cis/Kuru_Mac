      SUBROUTINE VIFLOW(VECTF,VECFL,ABETA,DMATP,HARDS,NSTRS,NTYPE,
     .                  PROPS,VECTG,VECGL,DGVEC,NSTR1,
     .                  YIELD,PREYS,VISCO,EXPON,DBETA,TLAMD,DLAMD,DTIME,
     .                  IVIFL,AZABA,AZABB,AZABC,ALFAP,NFORZ,TLAMB,
     .                  LARGE,IFREN,NDIME,XJACM,XJA3M,XJACI,XJA3I,
     .                  RECRY,SRECR,ERECR,ASREC,ALREC,YIELO,EFFPL,EFFP1,
     .                  SPTOT,TRMTX,SIGMA,DSTRA,DSTRL)
C***********************************************************************
C
C**** THIS SUBROUTINE EVALUATES THE VISCOPLASTIC PARAMETER LAMBDA
C
C     TLAMD=LAMBDA PARAMETER
C     DLAMD=INCREMENT OF LAMBDA PARAMETER
C     TLAMB=RESIDUAL OF LAMBDA PARAMETER (USEFUL TO CHECK CONVERGENCE)
C     ABETA=MULTIPLICATOR OF THE PSEUDO VISCOPLASTIC CONST. TENSOR
C     DBETA=not used !!
C
C     Notes:
C           Viscoelastic model is only possible with IVIFL=1, 4 & 6
C           The use of alfap ne 1 should be revised
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION VECTF(6), DMATP(6,6), DGVEC(6), PROPS(*),
     .          VECTG(6), VECFL(*),   VECGL(*), XJACM(NDIME,*), RECRY(*)
      DIMENSION SPTOT(*)
      DIMENSION TRMTX(NDIME,*), VECGX(6), VECFX(6),
     .          SIGMA(*), DSTRA(*), DSTRL(*), DSTRX(6), XJACI(NDIME,*)
C
C**** ROTATES (BACK) PLASTIC FLUX & NORMAL TENSORS FROM LOCAL TO GLOBAL
C     (MATERIAL) SYSTEM OF COORDINATES (THIS IS NEEDED FOR ORTHOTROPIC
C     MATERIALS)
C
      DO ISTR1=1,NSTR1
       VECGX(ISTR1)=VECGL(ISTR1)
       VECFX(ISTR1)=VECFL(ISTR1)
      ENDDO
      ISOTR=INT(PROPS(1))   ! 0=isotropic; 1=orthotropic
      IF(ISOTR.NE.0) THEN
       IF(IFREN.EQ.2.OR.IFREN.EQ.3) THEN
        XJA3M=1.0D0
        DETJX=1.0D0
        CALL PUFOBA(VECGL,VECGX,TRMTX,XJA3M,DETJX,    2)
        CALL PUFOBA(VECFL,VECFX,TRMTX,XJA3M,DETJX,    2)
       ELSE
        CALL RUNEND('ERROR: IFREN NE 2 OR 3 NOT IMPLEMENTED')
       ENDIF
      ENDIF
C
      IF(LARGE.EQ.0) THEN
       DO ISTR1=1,NSTR1
        VECTG(ISTR1)=VECGX(ISTR1)
        VECTF(ISTR1)=VECFX(ISTR1)
       ENDDO
      ELSE
C
C**** TRANSFORMS BACK PLASTIC FLUX & NORMAL TENSORS
C
       IF(IFREN.EQ.1.OR.IFREN.EQ.4) THEN                     ! no change
        DETJX=1.0D0
        CALL PUFOBA(VECGX,VECTG,XJACM,XJA3M,DETJX,    3)
        CALL PUFOBA(VECFX,VECTF,XJACM,XJA3M,DETJX,    3)
       ENDIF
       IF(IFREN.EQ.2.OR.IFREN.EQ.3.OR.IFREN.EQ.5.OR.
     .    IFREN.EQ.6.OR.IFREN.EQ.7.OR.IFREN.EQ.8.OR.
     .    IFREN.EQ.9.OR.IFREN.EQ.10) THEN
        DETJX=1.0D0
        CALL PUFOBA(VECGX,VECTG,XJACM,XJA3M,DETJX,    2)
        CALL PUFOBA(VECFX,VECTF,XJACM,XJA3M,DETJX,    2)
       ENDIF
      ENDIF        ! large.eq.0
C
      IBETA=0      ! initialization
      DLAMD=0.0D0
      ABETA=0.0D0
      TLAMB=0.0D0
C
C**** COMPUTES "DLAMD" ACCORDING WITH THE FLOW MODEL (VISCOUS LAW)
C
      IF(IVIFL.EQ.1) THEN
       IF(EXPON.EQ.1.0D0) THEN
        IF(YIELD.GT.PREYS) THEN
         DNUME=(YIELD-PREYS)/DTIME-VISCO*TLAMD/DTIME
         HARVI=VISCO/DTIME
         TLAMB=(YIELD-PREYS)-VISCO*TLAMD
         IBETA=1
        ENDIF
       ELSE
        IF(YIELD.GT.PREYS) THEN
         DAUX1=VISCO**EXPON
         DAUX2=DTIME*EXPON*(YIELD-PREYS)**(EXPON-1.0D0)
         DAUX3=(YIELD-PREYS)/(DTIME*EXPON)
         DNUME=DAUX3-DAUX1/DAUX2*TLAMD
         IF(VISCO.EQ.0.0D0.OR.EXPON.LT.1.0D0)     ! relaxation algorithm
     .    DNUME=DNUME*EXPON
         HARVI=DAUX1/DAUX2
         TLAMB=(YIELD-PREYS)-(VISCO**EXPON)*TLAMD/
     .                                    ((YIELD-PREYS)**(EXPON-1.0D0))
         IBETA=1
        ENDIF
       ENDIF
      ENDIF
C
      IF(IVIFL.EQ.2) THEN     ! to be revised
       IF(PREYS.GT.0.0D0) THEN
        ABETA=(EXP(EXPON*(YIELD-PREYS)/PREYS)-1.0D0)/VISCO
        DBETA=EXPON/PREYS*EXP(EXPON*(YIELD-PREYS)/PREYS)/VISCO
       ELSE
        CALL RUNEND('ERROR=TOT. HARDENING FUNC. LT 0.0')
       ENDIF
      ENDIF
C
      IF(IVIFL.EQ.3) THEN     ! to be revised
       IF(PREYS.GT.0.0D0) THEN
        ABETA=(((YIELD-PREYS)/PREYS)**EXPON)/VISCO
        DBETA=EXPON/PREYS*((YIELD-PREYS)/PREYS)*(EXPON-1.0D0)
        CALL RUNEND('ERROR=TOT. HARDENING FUNC. LT 0.0')
       ENDIF
      ENDIF
C
C**** ZABARAS' LAW
C
C     Idem IVIFL=1 with other definition of YIELD, PREYS & VISCO
C
C     YIELD=YIEL1=sinh(B sigma-bar)
C     PREYS=PREY1=0.0
C     VISCO=(1/(A e-C/T))**(1/n)
C
C     => lambda=(YIEL1/VISCO)**EXPON
C
C     Two forms: TANGENT (NFORZ=1) & SECANT (NFORZ=2)
C
C     The secant form is obatained by assuming:
C     partial (sinh(b sigma-bar))**n / partial sigma-bar = 
C     n (sinh(b sigma-bar))**(n-1)
C
      IF(IVIFL.EQ.4) THEN      ! to be revised (TLAMB?)
       IF(EXPON.EQ.0.0D0)
     .  CALL RUNEND('ERROR: EXPON=0 IN VISCOUS LAW')
       VISCO=AZABA*EXP(-AZABC)
       VISCO=(1.0D0/VISCO)**(1.0D0/EXPON)
       YIEL1=DSINH(AZABB*YIELD)
c      PREYS=0.0D0             ! to be used in the convergence criterion
       PREY1=0.0D0
       DEFAC=AZABB*DCOSH(AZABB*YIELD)
       IF(NFORZ.EQ.2) DEFAC=1.0D0
       IF(EXPON.EQ.1.0D0) THEN
        IF(YIEL1.GT.PREY1) THEN
         DNUME=(YIEL1-PREY1)/(DTIME*DEFAC)-VISCO*TLAMD/(DTIME*DEFAC)
         HARVI=VISCO/(DTIME*DEFAC)
        ENDIF
       ELSE
        IF(YIEL1.GT.PREY1) THEN
         DAUX1=VISCO**EXPON
         DAUX2=DTIME*EXPON*(YIEL1-PREY1)**(EXPON-1.0D0)*DEFAC
         IF(DAUX2.LT.0.0D0) DAUX2=-DAUX2
         DAUX3=(YIEL1-PREY1)/(DTIME*EXPON*DEFAC)
         DNUME=DAUX3-DAUX1/DAUX2*TLAMD
         HARVI=DAUX1/DAUX2
        ENDIF
       ENDIF
      ENDIF           ! ivifl.eq.4
C
C**** RECRYSTALLIZATION MODEL
C
      IF(IVIFL.EQ.5) THEN
       FT=RECRY(1)/RECRY(21)
       FL=(RECRY(8)/ERECR)**RECRY(2)
       DFLL=-FL*RECRY(2)/ERECR
       A1=DSINH(YIELD/(RECRY(5)*SRECR))
       E1=1.0D0/RECRY(6)
       FS=A1**E1
       A2=DCOSH(YIELD/(RECRY(5)*SRECR))
       E2=E1-1.0D0
       DFST=E1*A1**E2*A2/(RECRY(5)*SRECR)
       DFSS=-DFST*YIELD/SRECR
       BB=DTIME*FT*DFST
C
       DNUME=FT*FL*FS/BB-TLAMD/BB
       IF(DNUME.LT.0.0D0) DNUME=0.0D0  ! simplification to achieve conv.
       HARVI=1.0D0/BB-FT*DFLL*ALREC*DTIME/BB-FT*DFSS*ASREC*DTIME/BB
       TLAMB=DNUME*DTIME
       IBETA=1
      ENDIF           ! ivifl.eq.5
C
C*** SMA MODEL
C
      IF(IVIFL.EQ.6) THEN                        ! rate-independent law
       DYIEL=(YIELD-YIELO)/DTIME                 ! equiv. stress rate
       EFFPL=EFFPL-EFFP1                         ! correction
       IF(DYIEL.GT.0.0D0) THEN                   ! A to S transformation
        IF(YIELD.GE.VISCO.AND.YIELD.LT.AZABA) THEN
         DPSDS=-(1.0D0-EFFPL/EXPON)              ! >0
     .                             /(YIELD-AZABA)
         IF(YIELO.LT.VISCO)                      ! to avoid fict. return
     .    DYIEL=(YIELD-VISCO)/DTIME
         FUASA=DPSDS*DYIEL                       ! >0
         DNUME=(EXPON*FUASA-TLAMD)/(DPSDS*EXPON)
         HARVI=1.0D0/(DPSDS*EXPON)               ! >0
         TLAMB=EXPON*FUASA-TLAMD
         EFFPL=EFFPL+FUASA*DTIME*EXPON           ! correction
         IF(EFFPL.GT.EXPON) EFFPL=EXPON          ! upper bound
         IBETA=1
        ENDIF
        IF(YIELD.GE.AZABA.AND.EFFPL.LT.EXPON) THEN
         FXX=YIELD/YIELO                         ! guess
         DPSDS=-(1.0D0-EFFPL/EXPON)              ! >0 (factor estimated)
     .                             /(YIELD-AZABA*FXX)
         FUASA=DPSDS*DYIEL                       ! >0
         DNUME=(EXPON*FUASA-TLAMD)/(DPSDS*EXPON)
         HARVI=1.0D0/(DPSDS*EXPON)               ! >0
         TLAMB=EXPON*FUASA-TLAMD
         EFFPL=EXPON                             ! correction
         IBETA=1
        ENDIF
       ENDIF                                     ! dyiel.gt.0
C
       IF(DYIEL.LT.0.0D0) THEN                   ! S to A transformation
        IF(YIELD.LE.AZABC.AND.EFFPL.GT.0.0D0) THEN
         FXX=YIELD/YIELO                         ! guess
         DPSDS=EFFPL/EXPON/(YIELD-AZABC*FXX)     ! >0 (factor estimated)
         FUASA=DPSDS*DYIEL                       ! <0
         DNUME=(EXPON*FUASA-TLAMD)/(DPSDS*EXPON)
         HARVI=1.0D0/(DPSDS*EXPON)               ! >0
         TLAMB=EXPON*FUASA-TLAMD
         EFFPL=0.0D0                             ! correction
         DO ISTR1=1,NSTR1
          SPTOT(ISTR1)=0.0D0
         ENDDO
         IBETA=1
        ENDIF
        IF(YIELD.GT.AZABC.AND.YIELD.LE.AZABB) THEN
         DPSDS=EFFPL/EXPON/(YIELD-AZABC)         ! >0
         IF(YIELO.GT.AZABB)                      ! to avoid fict. return
     .    DYIEL=(YIELD-AZABB)/DTIME
         FUASA=DPSDS*DYIEL                       ! <0
         DNUME=(EXPON*FUASA-TLAMD)/(DPSDS*EXPON)
         HARVI=1.0D0/(DPSDS*EXPON)               ! >0
         TLAMB=EXPON*FUASA-TLAMD
         EFFPL=EFFPL+FUASA*DTIME*EXPON           ! correction
         IF(EFFPL.LT.0.0D0) THEN
          EFFPL=0.0D0                            ! lower bound
          DO ISTR1=1,NSTR1
           SPTOT(ISTR1)=0.0D0
          ENDDO
         ENDIF
         IBETA=1
        ENDIF
       ENDIF                                     ! dyiel.lt.0
      ENDIF           ! ivifl.eq.6
C
      IF(IVIFL.EQ.7) THEN                        ! Danilo's model
       IF(YIELD.GT.PREYS) THEN
c       EFFPX=EFFP1/DTIME    ! effective plastic deformation rate (epdr)
        IF(LARGE.EQ.0) THEN
         DO ISTR1=1,NSTR1
          DSTRL(ISTR1)=DSTRA(ISTR1)
         ENDDO
        ELSE
C
C**** TRANSFORMS INCREMENTAL DEFORMATION
C
         DETJX=1.0D0
         IF(IFREN.EQ.1.OR.IFREN.EQ.4) THEN                   ! no change
          CALL PUFOBA(DSTRA,DSTRL,XJACI,XJA3I,DETJX,    3)
         ENDIF
         IF(IFREN.EQ.2.OR.IFREN.EQ.3.OR.IFREN.EQ.5.OR.
     .      IFREN.EQ.6.OR.IFREN.EQ.7.OR.IFREN.EQ.8.OR.
     .      IFREN.EQ.9.OR.IFREN.EQ.10) THEN
          CALL PUFOBA(DSTRA,DSTRL,XJACI,XJA3I,DETJX,    2)   ! Almansi
         ENDIF
        ENDIF        ! large.eq.0
C
C**** ROTATES INCREMENTAL DEFORMATION FROM GLOBAL TO LOCAL (MATERIAL)
C     SYSTEM OF COORDINATES (THIS IS NEEDED FOR ORTHOTROPIC MATERIALS)
C
        IF(ISOTR.NE.0) THEN
         DO ISTR1=1,NSTR1
          DSTRX(ISTR1)=DSTRL(ISTR1)
         ENDDO
         IF(IFREN.EQ.2.OR.IFREN.EQ.3) THEN
          XJA3I=1.0D0
          DETJX=1.0D0
          CALL PUFOBA(DSTRX,DSTRL,TRMTX,XJA3I,DETJX,    1)
         ELSE
          CALL RUNEND('ERROR: IFREN NE 2 OR 3 NOT IMPLEMENTED')
         ENDIF
        ENDIF
C
        EFFPX=0.0D0          ! effective total deformation rate (etdr)
        DO ISTRS=1,NSTRS
         EFFPX=EFFPX+SIGMA(ISTRS)*DSTRL(ISTRS)/DTIME
        END DO
        EFFPX=DABS(EFFPX/YIELD)
C
        EFFPM=(AZABC/AZABA)**(-1.0D0/AZABB)      ! minimum etdr
        VISCO=AZABC
        IF(EFFPX.GT.EFFPM) VISCO=AZABA*DEXP(-AZABB*DLOG(EFFPX))
        DNUME=(YIELD-PREYS)/DTIME-VISCO*TLAMD/DTIME
        HARVI=VISCO/DTIME
        TLAMB=(YIELD-PREYS)-VISCO*TLAMD
        IBETA=1
       ENDIF
      ENDIF           ! ivifl.eq.7
C
C**** COMPUTE "DGVEC=(DMATP*VECTG)" PRODUCT
C
      IF(IBETA.EQ.1) THEN
       CALL PLPDXV(VECTG,DMATP,DGVEC,NSTRS,NTYPE,ISOTR,NSTR1)
C
C**** COMPUTE THE DENOMINATOR OF CONSISTENT PLASTIC FACTOR
C
       DENOM=ALFAP*HARDS+HARVI
       DO ISTRS=1,NSTRS
        DENOM=DENOM+ALFAP*VECTF(ISTRS)*DGVEC(ISTRS)
       ENDDO
       DLAMD=DNUME/DENOM
       ABETA=1.0D0/DENOM
      ENDIF           ! ibeta.eq.1
C
      RETURN
      END
