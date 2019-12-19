      SUBROUTINE LDSF01(ELCOD,LNODS,PROPS,THICK,
     .                  ALOAD,DERIV,ELEDG,GVECR,NOPRS,POSGP,PRESS,
     .                  PVECT,SHAPE,WEIGP,XJACM,ITASK,ATASK)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE EQUIVALENT NODAL FORCES DUE TO THE
C     TRACTION FORCES ACTING ON THE ELEMENT SURFACE
C     ( FOR ELEMENT NO. 1 )
C
C
C     Notes:
C
C     1) For triangular elements (NDIME=2), the equivalente nodal forces
C     are evaluated using the shape functions corresponding to the
C     quadrilateral elements considering the same kind of interpolation
C     functions (linear or parabolic) on the surface for both elements.
C     2) The local tangent axis (one in 2D and two in 3D) is defined by
C     the element numbering at the surface (boundary). The perpendicular
C     axis to the tangent gives the local normal axis (using, of course,
C     the right-hand convention). Note that the local axis are ordered
C     in the following form: normal & tangent.
C     3) The subroutine rulepw.f can not be called with a fixed NRULE
C     value because such value is reasigned inside it.
C     4) ITASK=0 => evaluates deformation-independent external load
C                   (ITASK=8 in elm030.f)
C        ITASK=1 => evaluates deformation-dependent external load
C                   (ITASK=4 in elm030.f)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
      INCLUDE 'inpo_om.f'
C
      COMMON/LDFILE/ITAPE
      COMMON/LDFILE1/IAXES
C
      DIMENSION ELCOD(NDIME,*),   LNODS(*),         PROPS(*)
      DIMENSION ALOAD(*),         DERIV(NDIME-1,*), ELEDG(NDIME,*),
     .          GVECR(*),         NOPRS(*),         POSGP(NDIME-1,*),
     .          PRESS(NNODE,*),   PVECT(*),         SHAPE(*), 
     .          WEIGP(*),         XJACM(NDIME,*)
C
      TWOPI=6.283185307179586D0
C
      IF(ITASK.EQ.1) THEN
       IF(NOPRS(1).EQ.0) RETURN  ! one node is enough to check (KLDSF=0)
      ENDIF
C
      IF(NDIME.EQ.3) GOTO 1000
C
C#######################################################################
C
C**** 2-D CASE
C
C#######################################################################
C
      NDIMB=1
      NGASS=2
      NODEG=3
      IF(NNODL.LE.4) NODEG=2
C
C**** IDENTIFY SAMPLING POINTS AND WEIGHTS
C
      CALL RULGAU(    2,POSGP,WEIGP)
C
C**** READ DATA LOCATING THE LOADED EDGE AND APPLIED LOAD
C
      IF(ITASK.EQ.0) THEN
       NPRIN=1
       CALL LISTEN('LDSF01',NPRIN,ITAPE)
       DO IODEG=1,NODEG
        NOPRS(IODEG)=INT(PARAM(IODEG))
        IF(NOPRS(IODEG).LE.0.OR.NOPRS(IODEG).GT.NPOIN)
     .   CALL RUNEND('LDSF01: WRONG NODE NUMBER')
       ENDDO
C
       WRITE(LURES,915) IELEM,IAXES,(NOPRS(IODEG),IODEG=1,NODEG)
C
       CALL LISTEN('LDSF01',NPRIN,ITAPE)
       ICOUN=0
       DO IODEG=1,NODEG
        DO IDIME=1,NDIME
         ICOUN=ICOUN+1
         PRESS(IODEG,IDIME)=PARAM(ICOUN)
        ENDDO
       ENDDO
C
       WRITE(LURES,920) ((PRESS(IODEG,IDIME),IDIME=1,NDIME),
     .                                       IODEG=1,NODEG)
      ENDIF              ! itask.eq.0
C
C**** STORE THE COORDINATES OF THE NODES OF THE ELEMENT EDGE
C
      DO IODEG=1,NODEG
       LNODE=NOPRS(IODEG)
       DO INODE=1,NNODL
        JNODE=LNODS(INODE)
        IF(LNODE.EQ.JNODE) THEN
         DO IDIME=1,NDIME
          ELEDG(IDIME,IODEG)=ELCOD(IDIME,INODE)
         ENDDO
        ENDIF
       ENDDO
      ENDDO
C
C**** ENTER LOOP FOR LINEAR NUMERICAL INTEGRATION
C
      DO 30 IGAUS=1,NGASS
C
      EXISP=POSGP(1,IGAUS)
c     IF(NNODL.EQ.3.OR.NNODL.EQ.6)
c    . EXISP=(EXISP+1.0D+00)/2.0D+00                      ! see note 1)
C
C**** EVALUATE THE SHAPE FUNCTIONS AT THE SAMPLING POINTS
C
      CALL SHAFUN(DERIV,EXISP,ETASP,EZETA,NDIMB,NODEG,NQUTR,    0,SHAPE)
C
C**** COMPUTE JACOBIAN
C
      DO IDIME=1,NDIME*NDIME
       XJACM(IDIME,1)=0.0D0
      ENDDO
C
      DO IDIME=1,NDIME
       DO IODEG=1,NODEG
        XJACM(1,IDIME)=XJACM(1,IDIME)+ELEDG(IDIME,IODEG)*
     .                 DERIV(1,IODEG)
       ENDDO
      ENDDO
C
      DAREA=XJACM(1,1)*XJACM(1,1)+XJACM(1,2)*XJACM(1,2)
      DAREA=DSQRT(DAREA)
C
C**** CALCULATE THE PRESSURE COMPONENTS AT THE GAUSSIAN POINTS 
C
      DO IDIME=1,NDIME
       GVECR(IDIME)=0.0D0
       DO IODEG=1,NODEG
        GVECR(IDIME)=GVECR(IDIME)+PRESS(IODEG,IDIME)*SHAPE(IODEG)
       ENDDO
      ENDDO
C
C**** COMPUTE LOCAL CARTESIAN SYSTEM, IF NECESSARY (IAXES=0)
C
      IF(IAXES.EQ.0) THEN
       XJA11=XJACM(1,1)
       XJA12=XJACM(1,2)
       XJACM(1,1)=-XJA12/DAREA   ! normalize; see note 2)
       XJACM(1,2)= XJA11/DAREA
       XJACM(2,1)= XJA11/DAREA
       XJACM(2,2)= XJA12/DAREA
      ELSE
       DO IDIME=1,NDIME
        DO JDIME=1,NDIME
         XJACM(IDIME,JDIME)=0.0D0
         IF(IDIME.EQ.JDIME) XJACM(IDIME,JDIME)=1.0D0
        ENDDO
       ENDDO
      ENDIF
C
C**** TRANSFORM PRESSURE COMPONENTS TO THE GLOBAL SYSTEM
C
      DO IDIME=1,NDIME
       PVECT(IDIME)=0.0D0
       DO JDIME=1,NDIME
        PVECT(IDIME)=PVECT(IDIME)+GVECR(JDIME)*XJACM(JDIME,IDIME)
       ENDDO
      ENDDO
c     PVECT(1)=GVECR(2)*XJACM(1,1)-GVECR(1)*XJACM(1,2)        ! old form
c     PVECT(2)=GVECR(1)*XJACM(1,1)+GVECR(2)*XJACM(1,2)
C
C**** COMPUTE ELEMENT AREA
C
      DAREA=DAREA*WEIGP(IGAUS)
c     IF(NNODL.EQ.3.OR.NNODL.EQ.6) DAREA=DAREA/2.0D+00     ! see note 1)
      IF(NTYPE.EQ.4) THICK=1.0D0               ! already set in inpset.f
      DAREA=DAREA*THICK
      IF(NTYPE.EQ.3) THEN
       RADUS=0.0D0
       DO IODEG=1,NODEG
        RADUS=RADUS+SHAPE(IODEG)*ELEDG(1,IODEG)
       ENDDO
       DAREA=DAREA*TWOPI*RADUS
      ENDIF
C
C**** COMPUTE ELEMENT CONTRIBUTION
C
      DO JODEG=1,NODEG
       LNODE=NOPRS(JODEG)
       DO INODE=1,NNODL
        KNODE=LNODS(INODE)
        IF(LNODE.EQ.KNODE) THEN
         SAREA=SHAPE(JODEG)*DAREA
         IEVAB=(INODE-1)*NDIME
         DO IDIME=1,NDIME
          IEVAB=IEVAB+1
          ALOAD(IEVAB)=ALOAD(IEVAB)+ATASK*PVECT(IDIME)*SAREA
         ENDDO
        ENDIF
       ENDDO
      ENDDO
C
   30 CONTINUE
C
      RETURN
C
C#######################################################################
C
C**** 3-D CASE
C
C#######################################################################
C
 1000 CONTINUE
C
      NDIMB=2
      IF(NNODL.EQ. 4) NODFC=3 
      IF(NNODL.EQ. 5) NODFC=3 
      IF(NNODL.EQ. 8) NODFC=4 
      IF(NNODL.EQ.10) NODFC=6 
      IF(NNODL.EQ.20) NODFC=8 
      IF(NNODL.EQ.27) NODFC=9
C
C**** IDENTIFY SAMPLING POINTS AND WEIGHTS
C
      IF(NNODL.EQ.4.OR.NNODL.EQ.5.OR.NNODL.EQ.10) THEN
       NGASS=3
       NRULX=3
       CALL RULEPW(    2,    3,NRULX,POSGP,WEIGP)          ! see note 3)
      ELSE
       NGASS=4
       NRULX=1
       CALL RULEPW(    2,    4,NRULX,POSGP,WEIGP)          ! see note 3)
      ENDIF
C
C**** READ DATA LOCATING THE LOADED FACE AND APPLIED LOAD
C
      IF(ITASK.EQ.0) THEN
       NPRIN=1
       CALL LISTEN('LDSF01',NPRIN,ITAPE)
       DO IODFC=1,NODFC
        NOPRS(IODFC)=INT(PARAM(IODFC))
        IF(NOPRS(IODFC).LE.0.OR.NOPRS(IODFC).GT.NPOIN)
     .   CALL RUNEND('LDSF01: WRONG NODE NUMBER')
       ENDDO
C
       WRITE(LURES,915) IELEM,IAXES,(NOPRS(IODFC),IODFC=1,NODFC)
C
       CALL LISTEN('LDSF01',NPRIN,ITAPE)
       ICOUN=0
       DO IODFC=1,NODFC
        DO IDIME=1,NDIME
         ICOUN=ICOUN+1
         PRESS(IODFC,IDIME)=PARAM(ICOUN)
        ENDDO
       ENDDO
C
       WRITE(LURES,920) ((PRESS(IODFC,IDIME),IDIME=1,NDIME),
     .                                       IODFC=1,NODFC)
      ENDIF              ! itask.eq.0
C
C**** STORE THE COORDINATES OF THE NODES OF THE ELEMENT FACE
C
      DO IODFC=1,NODFC
       LNODE=NOPRS(IODFC)
       DO INODE=1,NNODL
        JNODE=LNODS(INODE)
        IF(LNODE.EQ.JNODE) THEN
         DO IDIME=1,NDIME
          ELEDG(IDIME,IODFC)=ELCOD(IDIME,INODE)
         ENDDO
        ENDIF
       ENDDO
      ENDDO
C
C**** ENTER LOOP FOR LINEAR NUMERICAL INTEGRATION
C
      DO 130 IGAUS=1,NGASS
C
      EXISP=POSGP(1,IGAUS)
      ETASP=POSGP(2,IGAUS)
C
C**** EVALUATE THE SHAPE FUNCTIONS AT THE SAMPLING POINTS
C
      CALL SHAFUN(DERIV,EXISP,ETASP,EZETA,NDIMB,NODFC,NQUTR,    0,SHAPE)
C
C**** COMPUTE JACOBIAN
C
      DO IDIME=1,NDIME*NDIME
       XJACM(IDIME,1)=0.0D0
      ENDDO
C
      DO IDIME=1,NDIMB           ! compute two vectors on the surface
       DO JDIME=1,NDIME
        DO KNODE=1,NODFC
         XJACM(IDIME,JDIME)=XJACM(IDIME,JDIME)+
     .                      DERIV(IDIME,KNODE)*ELEDG(JDIME,KNODE)
        ENDDO
       ENDDO
      ENDDO
C
      VNOR1=0.0D0
      VNOR2=0.0D0
      VNOR3=0.0D0
      DO IDIMA=1,NDIME           ! compute normal as cross product
       IDIMB=IDIMA+1-IDIMA/3*3
       IDIMC=IDIMB+1-IDIMB/3*3
       XJACM(3,IDIMA)=XJACM(1,IDIMB)*XJACM(2,IDIMC)-
     .                XJACM(2,IDIMB)*XJACM(1,IDIMC)
       VNOR1=VNOR1+XJACM(1,IDIMA)*XJACM(1,IDIMA)
       VNOR3=VNOR3+XJACM(3,IDIMA)*XJACM(3,IDIMA)
      ENDDO
      DO IDIMA=1,NDIME           ! recompute 2nd vector as cross product
       IDIMB=IDIMA+1-IDIMA/3*3
       IDIMC=IDIMB+1-IDIMB/3*3
       XJACM(2,IDIMA)=XJACM(3,IDIMB)*XJACM(1,IDIMC)-
     .                XJACM(1,IDIMB)*XJACM(3,IDIMC)
       VNOR2=VNOR2+XJACM(2,IDIMA)*XJACM(2,IDIMA)
      ENDDO
      VNOR1=DSQRT(VNOR1)
      VNOR2=DSQRT(VNOR2)
      VNOR3=DSQRT(VNOR3)
      DAREA=VNOR3
C
C**** CALCULATE THE PRESSURE COMPONENTS AT THE GAUSSIAN POINTS 
C
      DO IDIME=1,NDIME
       GVECR(IDIME)=0.0D0
       DO IODFC=1,NODFC
        GVECR(IDIME)=GVECR(IDIME)+PRESS(IODFC,IDIME)*SHAPE(IODFC)
       ENDDO
      ENDDO
C
C**** COMPUTE LOCAL CARTESIAN SYSTEM, IF NECESSARY (IAXES=0)
C
      IF(IAXES.EQ.0) THEN
       DO IDIME=1,NDIME           ! normalize
        XJACM(1,IDIME)=XJACM(1,IDIME)/VNOR1
        XJACM(2,IDIME)=XJACM(2,IDIME)/VNOR2
        XJACM(3,IDIME)=XJACM(3,IDIME)/VNOR3
       ENDDO
       DO IDIME=1,NDIME
        XAUX1=XJACM(1,IDIME)
        XAUX2=XJACM(2,IDIME)
        XJACM(1,IDIME)=XJACM(3,IDIME)
        XJACM(2,IDIME)=XAUX1
        XJACM(3,IDIME)=XAUX2
       ENDDO
      ELSE
       DO IDIME=1,NDIME
        DO JDIME=1,NDIME
         XJACM(IDIME,JDIME)=0.0D0
         IF(IDIME.EQ.JDIME) XJACM(IDIME,JDIME)=1.0D0
        ENDDO
       ENDDO
      ENDIF
C
C**** TRANSFORM PRESSURE COMPONENTS TO THE GLOBAL SYSTEM 
C
      DO IDIME=1,NDIME
       PVECT(IDIME)=0.0D0
       DO JDIME=1,NDIME
        PVECT(IDIME)=PVECT(IDIME)+GVECR(JDIME)*XJACM(JDIME,IDIME)
       ENDDO
      ENDDO
C
C**** COMPUTE ELEMENT AREA
C
      DAREA=DAREA*WEIGP(IGAUS)
C
C**** COMPUTE ELEMENT CONTRIBUTION
C
      DO JODFC=1,NODFC
       LNODE=NOPRS(JODFC)
       DO INODE=1,NNODL
        KNODE=LNODS(INODE)
        IF(LNODE.EQ.KNODE) THEN
         SAREA=SHAPE(JODFC)*DAREA
         IEVAB=(INODE-1)*NDIME
         DO IDIME=1,NDIME
          IEVAB=IEVAB+1
          ALOAD(IEVAB)=ALOAD(IEVAB)+ATASK*PVECT(IDIME)*SAREA
         ENDDO
        ENDIF
       ENDDO
      ENDDO
C
  130 CONTINUE
C
      RETURN
c 900 FORMAT(10I5)
  915 FORMAT(1X,'ELEMENT NO. = ',I7,5X,'IAXES = ',I5,5X,'NODES =',8I7)
  920 FORMAT(6F10.3)
      END
