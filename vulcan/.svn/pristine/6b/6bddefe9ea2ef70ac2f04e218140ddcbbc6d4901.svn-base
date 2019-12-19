      SUBROUTINE MAMF32(DVOLU,PROPS,LNODS,SHAPE,WSTIF,EHIST,ELDIS,TENOD,
     .                  DISPL,ELCOD,VNORL,DMATX,AUXS1,AUXS2,BMATX,VTANL)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE STIFFNESS MATRIX FOR THE GAP ELEMENT
C     ( ELEMENT NO. 32 ) WHEN CONTACT ELEMENT SUBDIVISION IS CONSIDERED
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      COMMON/TUNING4/RITEN,RITEF,TRUPL,TRUPM,TOLGA,TOLGAM,PSUBC
C
      DIMENSION DVOLU(*),       PROPS(*),
     .          SHAPE(NNODN,*), WSTIF(*),
     .          EHIST(NHIST,*), ELDIS(NDOFC,*),
     .          TENOD(*),       DISPL(NDOFN,*),
     .          ELCOD(NDIME,*), VNORL(NDIME,*),
     .          LNODS(*)
C
      DIMENSION AUXS1(NEVAB,*), AUXS2(NEVAB,*),
     .          DMATX(NDIME,*), BMATX(NDIME,*), ! NDIME instead of NSTR1
     .          VTANL(NDIME,NDIME-1,*)
C
      NNOBO=NNODL/2
      IF(NOCOL.EQ.1) NNOBO=NNODL                   ! non-coincident mesh
      NEVBO=NNOBO*NDOFN
C
      RIGIN=PROPS(2)
      TEMPL=PROPS(5)
      INOTE=INT(PROPS(6))
      ICOMO=INT(PROPS(7))
      IITEN=INT(RITEN)      ! 0: tun. not considered; 1: tun. considered
      IITEF=INT(RITEF)
      NSUBC=INT(PSUBC)
C
C**** DEFINES JACOBIAN REGULARIZATION FACTORS & MINIMUM & MAXIMUM
C     ADDMISSIBLE CONTACT GAPS
C
      NNN20=20
      IF(ICONC.EQ.0) THEN
       IF(IITEN.EQ.0) THEN
        INP10=INT(PROPS(10))
        IF(INP10.EQ.1) THEN
         IITEF=10
         TRUPL=PROPS(11)
         TRUPM=PROPS(12)
         NNN20=INT(PROPS(13))
        ENDIF
        INP14=INT(PROPS(14))
        IF(INP14.EQ.1) THEN
         TOLGA=PROPS(15)
         TOLGAM=PROPS(16)
        ENDIF
        INP17=INT(PROPS(17))
        IF(INP17.EQ.1) THEN
         NSUBC=INT(PROPS(18))
        ENDIF
       ENDIF
      ENDIF
      IF(NSUBC.EQ.1) RETURN
C



      CALL MASS32(DVOLU,PROPS,LNODS,SHAPE,WSTIF,EHIST,ELDIS,TENOD,
     .            DISPL,ELCOD,VNORL,DMATX,AUXS1,AUXS2,BMATX,VTANL,   1)




C
      RETURN
      END
