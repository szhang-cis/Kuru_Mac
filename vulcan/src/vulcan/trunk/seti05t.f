      SUBROUTINE SETI05T(CARTDT,DVOLUT,ELCODT,EMASST,GPCODT,LNODST,
     .                   PROPST,SHAPET,THICKT,ELDIST,
     .                   CART1T,DVOL1T,GPCO1T,SHAP1T,
     .                   DERIVT,POSGPT,WEIGPT,XJACMT,
     .                   ELCO1T,EMAS1T,ELDI1T,DVOLIT,ELCOIT,DISPTT,
     .                   ITASIT,ITASKT)
C***********************************************************************
C
C**** THIS ROUTINE SETS UP SOME NEEDED CONSTANT MATRICES FOR
C     FUTURE USE ( ELEMENT NO. 5 )
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nued_om.f'   ! thermal-microstructural
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      COMMON/JACOBSTA/IERORT,KERORT
C
      DIMENSION CARTDT(NDIMET,NNODST,*), DVOLUT(*),
     .          ELCODT(NDIMET,*),        GPCODT(NDIMET,*),
     .          LNODST(*),               PROPST(*),
     .          SHAPET(NNODST,*),        ELDIST(NDOFCT,*),
     .          EMASST(*)
      DIMENSION DERIVT(NDIMET,NNODST,*), POSGPT(NDIMET,*),
     .          WEIGPT(*),               XJACMT(NDIMET,*),
     .          ELCO1T(NDIMET,*)
      DIMENSION CART1T(NDIMET,NNODST,*), DVOL1T(*),
     .          GPCO1T(NDIMET,*),        SHAP1T(NNODST,*),
     .          ELDI1T(NDOFCT,*),        EMAS1T(*)
      DIMENSION DVOLIT(*),               ELCOIT(NDIMET,*),
c    .          DISPTT(NDOFCM,*),        DISRTT(NDOFCM,*)
     .          DISPTT(NDOFCM,*)
C
C**** OBTAINS MATERIAL PROPERTIES
C
      CALL IDENPRT(PROPST)
      IF(IMICR.EQ.1) CALL IDEPROS(PROPST)
C
      IF(NMEMO1.EQ.0) THEN
       DO IDIMET=1,NDIMET
        DO INODLT=1,NNODLT
         ELCO1T(IDIMET,INODLT)=ELCODT(IDIMET,INODLT)
        ENDDO
       ENDDO
      ENDIF
C
      IF(NMEMO5.EQ.0) THEN
       DO IDOFCT=1,NDOFCT
        DO INODLT=1,NNODLT
         ELDI1T(IDOFCT,INODLT)=ELDIST(IDOFCT,INODLT)
        ENDDO
       ENDDO
      ENDIF
C
C**** DEALS WITH DEFORMED SHAPE
C
C     Notes:
C
C     LARGET=1 => ELCO1T: material coord. & ELCOIT=spatial coord.
C     LARGET=2 => ELCO1T: updated coord. & ELCOIT=spatial coord.
C     LARGET=3 => ELCO1T: spatial coord. & ELCOIT=material coord.
C
C     DISPTT: total displacements
C     DISRTT: incremental displacements
C
C     ITASIT=1 implies that seti05t.f is called more than once for the
C              same ITASKT
C
      IF(LARGET.NE.0) THEN
       IF(ITASIT.EQ.0) THEN
        IF(LARGET.EQ.1) THEN
         DO IDIMET=1,NDIMET
          DO INODLT=1,NNODLT
           ELCOIT(IDIMET,INODLT)=ELCO1T(IDIMET,INODLT)+
     .                           DISPTT(IDIMET,INODLT)
          ENDDO
         ENDDO
        ENDIF                    ! larget.eq.1
        IF(LARGET.EQ.2) THEN
         DO IDIMET=1,NDIMET
          DO INODLT=1,NNODLT
           ELCOIT(IDIMET,INODLT)=ELCO1T(IDIMET,INODLT)+
     .                           DISPTT(IDIMET,INODLT)
          ENDDO
         ENDDO
         DO IDIMET=1,NDIMET
          DO INODLT=1,NNODLT
           CALL RUNENDT('ERROR: LARGE=2 NOT IMPLEMENTED - SETI05T')
c          ELCO1T(IDIMET,INODLT)=ELCOIT(IDIMET,INODLT)-
c    .                           DISRTT(IDIMET,INODLT)
          ENDDO
         ENDDO
        ENDIF                    ! larget.eq.2
        IF(LARGET.EQ.3) THEN
         DO IDIMET=1,NDIMET
          DO INODLT=1,NNODLT
           ELCOIT(IDIMET,INODLT)=ELCO1T(IDIMET,INODLT)
          ENDDO
         ENDDO
         DO IDIMET=1,NDIMET
          DO INODLT=1,NNODLT
           ELCO1T(IDIMET,INODLT)=ELCOIT(IDIMET,INODLT)+
     .                           DISPTT(IDIMET,INODLT)
          ENDDO
         ENDDO
        ENDIF                   ! larget.eq.3
       ENDIF                    ! itasit.eq.0
      ENDIF                     ! larget.ne.0
C
      IF(NMEMO2.EQ.0) THEN
       IF(NRULET.GT.5)
     .  CALL RUNENDT('ERROR: NRULET GT 5 WITH NMEMO2=0')
       IF(ITASKT.EQ.1) THEN
        CALL SETM05T(CARTDT,DVOLUT,ELCO1T,EMASST,GPCODT,LNODST,
     .               PROPST,SHAPET,THICKT,
     .               DERIVT,POSGPT,WEIGPT,XJACMT,ITASKT)
       ELSE
        IF(LARGET.NE.0)                         ! just to compute DVOLIT
     .   CALL SETM05T(CART1T,DVOLIT,ELCOIT,EMAS1T,GPCO1T,LNODST,
     .                PROPST,SHAP1T,THICKT,
     .                DERIVT,POSGPT,WEIGPT,XJACMT,ITASKT)
        IF(NNODLT.EQ.NNODST) THEN
         DO IGAUST=1,NGAULT
          DVOL1T(IGAUST)=DVOLUT(IGAUST)
          DO INODLT=1,NNODLT
           SHAP1T(INODLT,IGAUST)=SHAPET(INODLT,IGAUST)
           DO IDIMET=1,NDIMET
            CART1T(IDIMET,INODLT,IGAUST)=CARTDT(IDIMET,INODLT,IGAUST)
            GPCO1T(IDIMET,IGAUST)=GPCODT(IDIMET,IGAUST)
           ENDDO
          ENDDO
         ENDDO
         DO INDEXT=1,NNODLT*NNODLT
          EMAS1T(INDEXT)=EMASST(INDEXT)
         ENDDO
        ELSE      ! shape, cartesian derivatives etc. must be recomputed
         CALL SETM05T(CART1T,DVOL1T,ELCO1T,EMAS1T,GPCO1T,LNODST,
     .                PROPST,SHAP1T,THICKT,
     .                DERIVT,POSGPT,WEIGPT,XJACMT,ITASKT)
        ENDIF                     ! nnodl.eq.nnods
        IF(LARGET.EQ.2.OR.LARGET.EQ.3)           ! recompute shape, etc.
     .   CALL SETM05T(CART1T,DVOL1T,ELCO1T,EMAS1T,GPCO1T,LNODST,
     .                PROPST,SHAP1T,THICKT,
     .                DERIVT,POSGPT,WEIGPT,XJACMT,ITASKT)
       ENDIF                      ! itaskt.eq.1
      ELSE
C
C**** NMEMO2=1 >> THERMAL PROBLEM INDEPENDENT OF ITASK
C
C     ITASK=1: checks the negative volumes
C     ITASK>1: compute all the necessary arrays
C
       IF(LARGET.NE.0)                          ! just to compute DVOLIT
     . CALL SETM05T(CART1T,DVOLIT,ELCOIT,EMAS1T,GPCO1T,LNODST,
     .              PROPST,SHAP1T,THICKT,
     .              DERIVT,POSGPT,WEIGPT,XJACMT,ITASKT)
       CALL SETM05T(CART1T,DVOL1T,ELCO1T,EMAS1T,GPCO1T,LNODST,
     .              PROPST,SHAP1T,THICKT,
     .              DERIVT,POSGPT,WEIGPT,XJACMT,ITASKT)
      ENDIF                       ! nmemo2.eq.0
C
      RETURN
      END
