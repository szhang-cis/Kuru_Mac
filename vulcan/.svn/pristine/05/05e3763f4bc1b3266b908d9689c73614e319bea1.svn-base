      SUBROUTINE SETI30(CARTD,DVOLU,ELCOD,EMASS,EPMTX,GPCOD,LNODS,
     .                  PROPS,RMAT1,SHAPE,STIFH,THICK,ELDIS,
     .                  CART1,DVOL1,GPCO1,SHAP1,BSBA1,
     .                  DERIV,POSGP,WEIGP,XJACM,RMAT2,CMEAN,SHAPR,
     .                  WSTIR,VNORL,ELCO1,EMAS1,BSBAR,ELDI1,
     .                  CARTS,SHAPS,GPCOS,
     .                  ELCO2,DISI2,DISP2,
     .                  STREI,VARII,
     .                  STRS0,EHIST,
     .                  ITASK)
C***********************************************************************
C
C**** THIS ROUTINE SETS UP SOME NEEDED CONSTANT MATRICES FOR
C     FUTURE USE ( ELEMENT NO. 30 )
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nued_om.f'   ! thermal-microstructural
      INCLUDE 'nuee_om.f'   ! mechanical-microstructural
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
      INCLUDE 'inpo_om.f'
C
      COMMON/TUNING30/PPART,PPARB,PPARI,ESTAB
      COMMON/JACOBSA/IEROR,KEROR
C
      DIMENSION CARTD(NDIME,NNODS,*), DVOLU(*),
     .          ELCOD(NDIME,*),       EPMTX(*),
     .          GPCOD(NDIME,*),       LNODS(*),
     .          PROPS(*),             RMAT1(NDIME,*), 
     .          SHAPE(NNODS,*),       EMASS(*),
     .          STIFH(*)
      DIMENSION DERIV(NDIME,NNODS,*), 
     .          POSGP(NDIME,*),       WEIGP(*),
     .          XJACM(NDIME,*),       RMAT2(NSTRS,*),
     .          ELCO1(NDIME,*)
      DIMENSION CMEAN(3,*),           SHAPR(NNODS,*),
     .          WSTIR(NNODS,*),       VNORL(NDIME,*),
     .          BSBAR(NSTR1,NEVAB,*)
      DIMENSION CART1(NDIME,NNODS,*), DVOL1(*),
     .          GPCO1(NDIME,*),       SHAP1(NNODS,*),
     .          BSBA1(NSTR1,NEVAB,*), EMAS1(*)
      DIMENSION ELDIS(NDOFC,*),       ELDI1(NDOFC,*)
      DIMENSION CARTS(NDIME,NNODS,*), SHAPS(NNODS,*), GPCOS(NDIME,*)
c     DIMENSION DVOLI(*),             ELCOI(NDIME,*)
      DIMENSION ELCO2(NDIME,*),
     .          DISI2(NDOFC,*),       DISP2(NDOFC,*)
      DIMENSION STREI(NPRE2,*),       VARII(NPRE3,*)
      DIMENSION STRS0(NSTR1,*),       EHIST(NHIST,*)
C
      TWOPI=6.283185307179586D0
C
C**** OBTAINS MATERIAL PROPERTIES
C
      CALL IDENPR(PROPS)
      IF(IMICRM.GT.0.OR.IMICR.GT.0) CALL IDEPROM(PROPS)
C
      IF(NMEMO1M.EQ.0) THEN
       DO IDIME=1,NDIME
        DO INODE=1,NNODL
         ELCO1(IDIME,INODE)=ELCOD(IDIME,INODE)
        ENDDO
       ENDDO
      ENDIF
C
      IF(NMEMO5M.EQ.0) THEN
       DO IDOFC=1,NDOFC
        DO INODE=1,NNODL
         ELDI1(IDOFC,INODE)=ELDIS(IDOFC,INODE)
        ENDDO
       ENDDO
      ENDIF
C
C**** CONTROL DEFORMED SHAPE (FOR SMALL & LARGE STRAINS)
C
      IF(ITASK.EQ.3.OR.ITASK.EQ.4) THEN    ! stiffness m. or residual f.
       IF(KPROB.NE.5) THEN                 ! exclude incompressib. cond.
        DO IDIME=1,NDIME
         DO INODE=1,NNODL
          ELCO1(IDIME,INODE)=ELCO1(IDIME,INODE)+ELDI1(IDIME,INODE)
         ENDDO
        ENDDO
C
        CALL SETM30(CART1,DVOL1,ELCO1,EMAS1,EPMTX,GPCO1,LNODS,
     .              PROPS,RMAT1,SHAP1,STIFH,THICK,
     .              DERIV,POSGP,WEIGP,XJACM,RMAT2,CMEAN,SHAPR,
     .              WSTIR,VNORL,BSBA1,CARTS,SHAPS,GPCOS,    1,ITASK)
C
        DO IDIME=1,NDIME
         DO INODE=1,NNODL
          ELCO1(IDIME,INODE)=ELCO1(IDIME,INODE)-ELDI1(IDIME,INODE)
         END DO
        ENDDO
C
        IF(IEROR.NE.0) RETURN
       ENDIF                             ! kprob.ne.5
      ENDIF                              ! itask.eq.3.or.itask.eq.4
C
C**** DEALS WITH DEFORMED SHAPE
C
C     Notes:
C
C     LARGE=1 => ELCO1: material coord. & ELCOI=spatial coord.
C     LARGE=2 => ELCO1: updated coord. & ELCOI=spatial coord.
C     LARGE=3 => ELCO1: spatial coord. & ELCOI=material coord.
C
C     DISPT: total displacements
C     DISRT: incremental displacements
C
c     IF(LARGE.NE.0) THEN
c      IF(LARGE.EQ.1) THEN
c       DO IDIME=1,NDIME
c        DO INODE=1,NNODL
c         ELCOI(IDIME,INODE)=ELCO1(IDIME,INODE)+
c    .                       DISPT(IDIME,INODE)
c        ENDDO
c       ENDDO
c      ENDIF                    ! large.eq.1
c      IF(LARGE.EQ.2) THEN
c       DO IDIME=1,NDIME
c        DO INODE=1,NNODL
c         ELCOI(IDIME,INODE)=ELCO1(IDIME,INODE)+
c    .                       DISPT(IDIME,INODE)
c        ENDDO
c       ENDDO
c       DO IDIME=1,NDIME
c        DO INODE=1,NNODL
c         CALL RUNEND('ERROR: LARGE=2 NOT IMPLEMENTED - SETM30')
c         ELCO1(IDIME,INODE)=ELCOI(IDIME,INODE)-
c    .                       DISRT(IDIME,INODE)
c        ENDDO
c       ENDDO
c      ENDIF                    ! large.eq.2
c      IF(LARGE.EQ.3) THEN
c       DO IDIME=1,NDIME
c        DO INODE=1,NNODL
c         ELCOI(IDIME,INODE)=ELCO1(IDIME,INODE)
c        ENDDO
c       ENDDO
c       DO IDIME=1,NDIME
c        DO INODE=1,NNODL
c         ELCO1(IDIME,INODE)=ELCOI(IDIME,INODE)+
c    .                       DISPT(IDIME,INODE)
c        ENDDO
c       ENDDO
c      ENDIF                   ! large.eq.3
c     ENDIF                    ! large.ne.0
C
C**** UPDATES GEOMETRY FOR LARGE DISPLACEMENTS
C
      IF(LARGE.NE.0) THEN
       DO IDIME=1,NDIME
        DO INODE=1,NNODL
         ELDI1(IDIME,INODE)=ELCO1(IDIME,INODE)+ELDI1(IDIME,INODE)
        END DO
       ENDDO
      ENDIF
C
      IF(NMEMO2M.EQ.0) THEN
       IF(ITASK.EQ.1) THEN
        CALL SETM30(CARTD,DVOLU,ELCO1,EMASS,EPMTX,GPCOD,LNODS,
     .              PROPS,RMAT1,SHAPE,STIFH,THICK,
     .              DERIV,POSGP,WEIGP,XJACM,RMAT2,CMEAN,SHAPR,
     .              WSTIR,VNORL,BSBAR,CARTS,SHAPS,GPCOS,    0,ITASK)
       ELSE
c       IF(LARGE.NE.0)                           ! just to compute DVOLI
c    .   CALL SETM30(CART1,DVOLI,ELCOI,EMAS1,EPMTX,GPCO1,LNODS,
c    .               PROPS,RMAT1,SHAP1,STIFH,THICK,
c    .               DERIV,POSGP,WEIGP,XJACM,RMAT2,CMEAN,SHAPR,
c    .               WSTIR,VNORL,BSBA1,CARTS,SHAPS,GPCOS,    0,ITASK)
        IF(NNODL.EQ.NNODS) THEN
         DO IGAUS=1,NGAUL
          DVOL1(IGAUS)=DVOLU(IGAUS)
          DO INODE=1,NNODL
           SHAP1(INODE,IGAUS)=SHAPE(INODE,IGAUS)
           DO IDIME=1,NDIME
            CART1(IDIME,INODE,IGAUS)=CARTD(IDIME,INODE,IGAUS)
            GPCO1(IDIME,IGAUS)=GPCOD(IDIME,IGAUS)
           ENDDO
          ENDDO
          DO ISTR1=1,NSTR1
           DO IEVAB=1,NEVAB
            BSBA1(ISTR1,IEVAB,IGAUS)=BSBAR(ISTR1,IEVAB,IGAUS)
           ENDDO
          ENDDO
         ENDDO
         DO INDEX=1,NNODL*NNODL 
          EMAS1(INDEX)=EMASS(INDEX)
         ENDDO
        ELSE      ! shape, cartesian derivatives etc. must be recomputed
         CALL SETM30(CART1,DVOL1,ELCO1,EMAS1,EPMTX,GPCO1,LNODS,
     .               PROPS,RMAT1,SHAP1,STIFH,THICK,
     .               DERIV,POSGP,WEIGP,XJACM,RMAT2,CMEAN,SHAPR,
     .               WSTIR,VNORL,BSBA1,CARTS,SHAPS,GPCOS,    0,ITASK)
        ENDIF              ! nnodl.eq.nnods
c       IF(LARGE.EQ.2.OR.LARGE.EQ.3)             ! recompute shape, etc.
c    .   CALL SETM30(CART1,DVOL1,ELCO1,EMAS1,EPMTX,GPCO1,LNODS,
c    .               PROPS,RMAT1,SHAP1,STIFH,THICK,
c    .               DERIV,POSGP,WEIGP,XJACM,RMAT2,CMEAN,SHAPR,
c    .               WSTIR,VNORL,BSBA1,CARTS,SHAPS,GPCOS,    0,ITASK)
       ENDIF               ! itask.eq.1
      ELSE
C
C**** NMEMO2M=1 >> MECHANICAL PROBLEM INDEPENDENT OF ITASK
C
C     ITASK=1: compute the normal vector
C     ITASK>1: compute all the necessary arrays
C
c      IF(LARGE.NE.0)                            ! just to compute DVOLI
c    . CALL SETM30(CART1,DVOLI,ELCOI,EMAS1,EPMTX,GPCO1,LNODS,
c    .             PROPS,RMAT1,SHAP1,STIFH,THICK,
c    .             DERIV,POSGP,WEIGP,XJACM,RMAT2,CMEAN,SHAPR,
c    .             WSTIR,VNORL,BSBA1,CARTS,SHAPS,GPCOS,    0,ITASK)
       CALL SETM30(CART1,DVOL1,ELCO1,EMAS1,EPMTX,GPCO1,LNODS,
     .             PROPS,RMAT1,SHAP1,STIFH,THICK,
     .             DERIV,POSGP,WEIGP,XJACM,RMAT2,CMEAN,SHAPR,
     .             WSTIR,VNORL,BSBA1,CARTS,SHAPS,GPCOS,    0,ITASK)
      ENDIF                ! nmemo2m.eq.0
C
C**** UPDATES GEOMETRY FOR DEFORMATION-DEPENDENT LOADS (LARGE > 0 IN
C     THIS CASE)
C
      IF(NLDSF.EQ.1) THEN
       IF(ILDSF.EQ.1) THEN
c       IF(ITASK.EQ.3.OR.ITASK.EQ.4) THEN
        IF(ITASK.EQ.4) THEN                  ! assumption (see proinp.f)
         DO IDIME=1,NDIME
          DO INODE=1,NNODL
           ELCO2(IDIME,INODE)=ELDI1(IDIME,INODE)
           IF(LLDSF.EQ.0) THEN               ! linearized computation
            IF(ITASK.EQ.4)              ! geometry is the same for R & K
     .       ELCO2(IDIME,INODE)=ELDI1(IDIME,INODE)-DISI2(IDIME,INODE)
            ELCO2(IDIME,INODE)=ELCO2(IDIME,INODE)-DISP2(IDIME,INODE)
           ENDIF
          END DO
         END DO
        ENDIF              ! itask.eq.4
       ENDIF               ! ildsf.eq.1
      ENDIF                ! nldsf.eq.1
C
      IF(INITV.EQ.1.AND.IPRCO.EQ.0) THEN    ! non-standard initial cond.
       IF(ITASK.EQ.1) THEN
        DO IGAUS=1,NGAUL
         IF(NPRE2.GT.0) THEN
          DO IPRE2=1,NPRE2                  ! initial stresses
           STREX=0.0D0
           DO INODE=1,NNODL
            STREX=STREX+SHAP1(INODE,IGAUS)*STREI(IPRE2,INODE)
           ENDDO
           STRS0(IPRE2,IGAUS)=STREX
          ENDDO
         ENDIF
         IF(NPRE3.GT.0) THEN
          DO IPRE3=1,NPRE3                  ! initial internal variables
           VARIX=0.0D0
           DO INODE=1,NNODL
            VARIX=VARIX+SHAP1(INODE,IGAUS)*VARII(IPRE3,INODE)
           ENDDO
           EHIST(IPLAS(3)-1+IPRE3,IGAUS)=VARIX     ! to be generalized !
          ENDDO
         ENDIF
        ENDDO              ! igaus=1,ngaul
       ENDIF               ! itask.eq.1
      ENDIF                ! initv.eq.1.and.iprco.eq.0
C
      RETURN
      END
