      SUBROUTINE SETDMG(DESIG,DMATX,DVOLU,ELCOD,NDIME,NNODE,
     .                  NSTR1,PRESG,PROPS,SGTOT,SIGMA,STRAN,RADUS,
     .                  VISCO,EQUIV,DAMAG,THRES,ALENG,AVALU,TDAMA,THICK)
C***********************************************************************
C
C****THIS ROUTINE DETERMINES IF THE POINT IS DAMAGED
C
C    Input variables:
C      
C      DMATX(NSTRS,NSTRS) Elastic matrix
C      DVOLU(NGAUL)       Volume of the integration points
C      ELCOD(NDIME,NNODE) Nodal coordinates
C      PRESG(NSTR1)       Previous total stresses
C      PROPS(NPROP)       Material properties 
C      SIGMA(NSTR1)       Total elastic stresses
C      STRAN(NSTR1)       Current total strains
C      RADUS              Radius ( axisymmetric )
C
C    Output variables:
C
C      TDAMA              Time for onset of damage
C
C    Input/Output variables:
C
C      ALENG(1:NDIME)     Characteristic lengths
C      AVALU              Softening parameter
C      EQUIV              Equivalent strain
C      DAMAG              Damage internal value
C      THRES              Threshold for damage
C      SGTOT(NSTR1)       Current total stresses
C      VISCO(1:NSTR1)     Internal variables for viscoelasticity
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION ALENG(*),           DESIG(*),  DMATX(*), 
     .          ELCOD(NDIME,*),     PRESG(*),  PROPS(*),  SGTOT(*), 
     .          SIGMA(*),           STRAN(*),  VISCO(*),  DVOLU(*)
      DIMENSION BETAM(9), STPRI(3)
      PARAMETER (TWOPI=6.283185307179586, ZEROM=1.0D-10)
C
C***DEFINE MATERIAL PROPERTIES
C
      YOUNG=PROPS( 2)
      FTULT=PROPS(32)
      GFVAL=PROPS(33)
      THRE0=FTULT/DSQRT(YOUNG)              ! initial threshold value
      IF(THRES.EQ.0.0D00) THRES=THRE0
C
C***DEAL WITH VISCOELASTICITY, IF NECESSARY
C
      RATIO=(YOUNG-PROPS(41))/YOUNG
      TIMER=PROPS(42)
      IF(TIMER.NE.0.0) THEN
C
        RATTI=DTIME/TIMER
        EXPON=EXP(-RATTI)
        UNOEX=1.0D00-EXPON
        COEFA=UNOEX*RATIO
C
        DO ISTRE=1,NSTR1
          PRESG(ISTRE)=((PRESG(ISTRE)/(1.0D00-DAMAG))+VISCO(ISTRE))*
     .                 RATIO
          VISCO(ISTRE)=COEFA*PRESG(ISTRE)+EXPON*VISCO(ISTRE)
          SIGMA(ISTRE)=SIGMA(ISTRE)-VISCO(ISTRE)
        ENDDO
C
      ENDIF
C
C***COMPUTE EQUIVALENT DAMAGE VARIABLE
C
      IF(NCRIT.EQ.21) THEN                    ! ISOTROPIC DAMAGE
C
        EQUIV=0.0D00
        DO ISTRS=1,NSTR1
          EQUIV=EQUIV+SIGMA(ISTRS)*SIGMA(ISTRS)
        ENDDO
        EQUIV=DSQRT(EQUIV/YOUNG)
C
      ELSE IF(NCRIT.EQ.22) THEN         ! TENSILE DAMAGE (stress based)
C
        CALL PRIVAL(NSTR1,SIGMA,STPRI)
C
        EQUIV=0.0D00
        DO INDEX=1,NDIME
          IF(STPRI(INDEX).GT.0.0D00) EQUIV=EQUIV+STPRI(INDEX)*
     .                                     STPRI(INDEX)
        ENDDO
        EQUIV=DSQRT(EQUIV/YOUNG)
C
      ELSE IF(NCRIT.EQ.23) THEN           ! MODEL III (AUSTRIA PAPER)
C
        CALL PRIVAL(NSTR1,SIGMA,STPRI)
C
        XNUME=0.0D00
        XDENO=0.0D00
        EQUIV=0.0D00
        DO INDEX=1,NDIME
          XDENO=XDENO+ABS(STPRI(INDEX))
          IF(STPRI(INDEX).GT.0.0D00) XNUME=XNUME+STPRI(INDEX)
          EQUIV=EQUIV+STPRI(INDEX)*STPRI(INDEX)
        ENDDO
        IF(XDENO.GT.1.0D-15) THEN
          THETA=XNUME/XDENO
          XMULT=THETA+(1.0D00-THETA)/10.0D00
          EQUIV=XMULT*DSQRT(EQUIV/YOUNG)
        ENDIF
C
      ENDIF
C
C***CHECK FOR DAMAGE & UPDATE DAMAGE INTERNAL VARIABLE 
C
      IF(EQUIV.GT.THRES) THEN             
C
C     (a) If necessary define initial parameters
C
      IF(ALENG(1).EQ.0.0D00) THEN
        TVOLU=0.0D00                             ! element volumen
        DO I=1,NGAUL
          TVOLU=TVOLU+DVOLU(I)
        ENDDO
        IF(NTYPE.NE.4) TVOLU=TVOLU/THICK
        IF(NTYPE.EQ.3) TVOLU=TVOLU/(TWOPI*RADUS*THICK)   ! not exact !!
C
        ALENG(1)=TVOLU**(1.0D00/DFLOAT(NDIME))
        DLENG=ALENG(1)
        GSPEC=GFVAL/DLENG                     ! softening parameter
        AINVE=(GSPEC/(THRE0*THRE0))-0.5D00
        IF(AINVE.LE.0.0D00) AINVE=1.0E-06
        AVALU=1.0D00/AINVE                   
      ENDIF
C
C    (b) update threshold for damage
C
        THRES=EQUIV
C
C    (c) update damage
C
        RATIO=THRES/THRE0
        EXPVA=AVALU*(1.0-RATIO)
        DAMAG=1.0D00-(EXP(EXPVA))/RATIO
        IF(DAMAG.EQ.1.0D00) DAMAG=1.0D00-ZEROM
C
C    (d) record time for onset of damage 
C
        IF(TDAMA.EQ.0.0D+00) TDAMA=TTIME
C
      ENDIF
C
C***COMPUTE REAL STRESSES ACCORDING TO VALUE OF DAMAGE
C
      REDUC=1.0D00-DAMAG
      DO ISTRE=1,NSTR1
        SGTOT(ISTRE)=REDUC*SIGMA(ISTRE)
      ENDDO
C
      RETURN
      END
