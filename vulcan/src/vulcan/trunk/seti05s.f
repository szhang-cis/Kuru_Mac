      SUBROUTINE SETI05S(CARTDS,DVOLUS,ELCODS,EMASSS,GPCODS,LNODSS,
     .                   PROPSS,SHAPES,THICKS,ELDISS,
     .                   CART1S,DVOL1S,GPCO1S,SHAP1S,
     .                   DERIVS,POSGPS,WEIGPS,XJACMS,
     .                   ELCO1S,EMAS1S,ELDI1S,DVOLIS,ELCOIS,DISPTS,
     .                   ITASIS,ITASKS)
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
      INCLUDE 'addi_oms.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nued_om.f'   ! thermal-microstructural
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
C
      COMMON/JACOBSSA/IERORS,KERORS
C
      DIMENSION CARTDS(NDIMES,NNODSS,*), DVOLUS(*),
     .          ELCODS(NDIMES,*),        GPCODS(NDIMES,*),
     .          LNODSS(*),               PROPSS(*),
     .          SHAPES(NNODSS,*),        ELDISS(NDOFCS,*),
     .          EMASSS(*)
      DIMENSION DERIVS(NDIMES,NNODSS,*), POSGPS(NDIMES,*),
     .          WEIGPS(*),               XJACMS(NDIMES,*),
     .          ELCO1S(NDIMES,*)
      DIMENSION CART1S(NDIMES,NNODSS,*), DVOL1S(*),
     .          GPCO1S(NDIMES,*),        SHAP1S(NNODSS,*),
     .          ELDI1S(NDOFCS,*),        EMAS1S(*)
      DIMENSION DVOLIS(*),               ELCOIS(NDIMES,*),
c    .          DISPTS(NDOFCMS,*),       DISRTS(NDOFCMS,*)
     .          DISPTS(NDOFCMS,*)
C
C**** OBTAINS MATERIAL PROPERTIES
C
      CALL IDENPRS(PROPSS)
c     IF(IMICR.EQ.1) CALL IDEPROS(PROPST)
C
      IF(NMEMO1S.EQ.0) THEN
       DO IDIMES=1,NDIMES
        DO INODLS=1,NNODLS
         ELCO1S(IDIMES,INODLS)=ELCODS(IDIMES,INODLS)
        ENDDO
       ENDDO
      ENDIF
C
      IF(NMEMO5S.EQ.0) THEN
       DO IDOFCS=1,NDOFCS
        DO INODLS=1,NNODLS
         ELDI1S(IDOFCS,INODLS)=ELDISS(IDOFCS,INODLS)
        ENDDO
       ENDDO
      ENDIF
C
C**** DEALS WITH DEFORMED SHAPE
C
C     Notes:
C
C     LARGES=1 => ELCO1S: material coord. & ELCOIS=spatial coord.
C     LARGES=2 => ELCO1S: updated coord. & ELCOIS=spatial coord.
C     LARGES=3 => ELCO1S: spatial coord. & ELCOIS=material coord.
C
C     DISPTS: total displacements
C     DISRTS: incremental displacements
C
C     ITASIS=1 implies that seti05t.f is called more than once for the
C              same ITASKT
C
      IF(LARGES.NE.0) THEN
       IF(ITASIS.EQ.0) THEN
        IF(LARGES.EQ.1) THEN
         DO IDIMES=1,NDIMES
          DO INODLS=1,NNODLS
           ELCOIS(IDIMES,INODLS)=ELCO1S(IDIMES,INODLS)+
     .                           DISPTS(IDIMES,INODLS)
          ENDDO
         ENDDO
        ENDIF                    ! larges.eq.1
        IF(LARGES.EQ.2) THEN
         DO IDIMES=1,NDIMES
          DO INODLS=1,NNODLS
           ELCOIS(IDIMES,INODLS)=ELCO1S(IDIMES,INODLS)+
     .                           DISPTS(IDIMES,INODLS)
          ENDDO
         ENDDO
         DO IDIMES=1,NDIMES
          DO INODLS=1,NNODLS
           CALL RUNENDS('ERROR: LARGE=2 NOT IMPLEMENTED - SETI05S')
c          ELCO1S(IDIMES,INODLS)=ELCOIS(IDIMES,INODLS)-
c    .                           DISRTS(IDIMES,INODLS)
          ENDDO
         ENDDO
        ENDIF                    ! larges.eq.2
        IF(LARGES.EQ.3) THEN
         DO IDIMES=1,NDIMES
          DO INODLS=1,NNODLS
           ELCOIS(IDIMES,INODLS)=ELCO1S(IDIMES,INODLS)
          ENDDO
         ENDDO
         DO IDIMES=1,NDIMES
          DO INODLS=1,NNODLS
           ELCO1S(IDIMES,INODLS)=ELCOIS(IDIMES,INODLS)+
     .                           DISPTS(IDIMES,INODLS)
          ENDDO
         ENDDO
        ENDIF                   ! larges.eq.3
       ENDIF                    ! itasis.eq.0
      ENDIF                     ! larges.ne.0
C
      IF(NMEMO2S.EQ.0) THEN
       IF(NRULES.GT.5)
     .  CALL RUNENDS('ERROR: NRULES GT 5 WITH NMEMO2S=0')
       IF(ITASKS.EQ.1) THEN
        CALL SETM05S(CARTDS,DVOLUS,ELCO1S,EMASSS,GPCODS,LNODSS,
     .               PROPSS,SHAPES,THICKS,
     .               DERIVS,POSGPS,WEIGPS,XJACMS)
       ELSE
        IF(LARGES.NE.0)                         ! just to compute DVOLIS
     .   CALL SETM05S(CART1S,DVOLIS,ELCOIS,EMAS1S,GPCO1S,LNODSS,
     .                PROPSS,SHAP1S,THICKS,
     .                DERIVS,POSGPS,WEIGPS,XJACMS)
        IF(NNODLS.EQ.NNODSS) THEN
         DO IGAUSS=1,NGAULS
          DVOL1S(IGAUSS)=DVOLUS(IGAUSS)
          DO INODLS=1,NNODLS
           SHAP1S(INODLS,IGAUSS)=SHAPES(INODLS,IGAUSS)
           DO IDIMES=1,NDIMES
            CART1S(IDIMES,INODLS,IGAUSS)=CARTDS(IDIMES,INODLS,IGAUSS)
            GPCO1S(IDIMES,IGAUSS)=GPCODS(IDIMES,IGAUSS)
           ENDDO
          ENDDO
         ENDDO
         DO INDEXS=1,NNODLS*NNODLS
          EMAS1S(INDEXS)=EMASSS(INDEXS)
         ENDDO
        ELSE      ! shape, cartesian derivatives etc. must be recomputed
         CALL SETM05S(CART1S,DVOL1S,ELCO1S,EMAS1S,GPCO1S,LNODSS,
     .                PROPSS,SHAP1S,THICKS,
     .                DERIVS,POSGPS,WEIGPS,XJACMS)
        ENDIF                     ! nnodl.eq.nnods
        IF(LARGES.EQ.2.OR.LARGES.EQ.3)           ! recompute shape, etc.
     .   CALL SETM05S(CART1S,DVOL1S,ELCO1S,EMAS1S,GPCO1S,LNODSS,
     .                PROPSS,SHAP1S,THICKS,
     .                DERIVS,POSGPS,WEIGPS,XJACMS)
       ENDIF                      ! itaskt.eq.1
      ELSE
C
C**** NMEMO2S=1 >> THERMAL PROBLEM INDEPENDENT OF ITASK
C
C     ITASK=1: checks the negative volumes
C     ITASK>1: compute all the necessary arrays
C
       IF(LARGES.NE.0)                          ! just to compute DVOLIS
     . CALL SETM05S(CART1S,DVOLIS,ELCOIS,EMAS1S,GPCO1S,LNODSS,
     .              PROPSS,SHAP1S,THICKS,
     .              DERIVS,POSGPS,WEIGPS,XJACMS)
       CALL SETM05S(CART1S,DVOL1S,ELCO1S,EMAS1S,GPCO1S,LNODSS,
     .              PROPSS,SHAP1S,THICKS,
     .              DERIVS,POSGPS,WEIGPS,XJACMS)
      ENDIF                       ! nmemo2s.eq.0
C
      RETURN
      END
