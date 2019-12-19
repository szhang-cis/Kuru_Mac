      SUBROUTINE SETI32(DVOLU,ELCOD,GPCOD,LNODS,PROPS,SHAPE,THICK,ELDIS,
     .                  VNORI,
     .                  DVOL1,GPCO1,SHAP1,
     .                  DERIV,POSGP,WEIGP,XJACM,
     .                  VNORL,ELCO1,ELDI1,STRSG,
     .                  DVOL2,ELCO2,GPCO2,SHAP2,DISI2,DISP2,VTANL,
     .                  ITASK)
C***********************************************************************
C
C**** THIS ROUTINE SETS UP SOME NEEDED CONSTANT MATRICES FOR
C     FUTURE USE ( ELEMENT NO. 32 )
C
C     Note: SHAPE, SHAP1, ELDIS, ELDI1 & DERIV must be dimensioned
C           with NNODN due to non-coincident contact mesh (NOCOL=1)
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
C
      COMMON/TUNING4/RITEN,RITEF,TRUPL,TRUPM,TOLGA,TOLGAM
C
      COMMON/JACOBSA/IEROR,KEROR
C
      DIMENSION DVOLU(*),             ELCOD(NDIME,*),
     .          GPCOD(NDIME,*),       LNODS(*),
     .          PROPS(*),             SHAPE(NNODN,*)
      DIMENSION DERIV(NDIME,NNODN,*), 
     .          POSGP(NDIME,*),       WEIGP(*),
     .          XJACM(NDIME,*),
     .          ELCO1(NDIME,*),       VNORI(NDIME,*)
      DIMENSION VNORL(NDIME,*)
      DIMENSION DVOL1(*),
     .          GPCO1(NDIME,*),       SHAP1(NNODN,*)
      DIMENSION ELDIS(NDOFC,*),       ELDI1(NDOFC,*)
      DIMENSION STRSG(NSTR1,*)          ! useful for non-coincident mesh
      DIMENSION DVOL2(*),             ELCO2(NDIME,*),
     .          GPCO2(NDIME,*),       SHAP2(NNODN,*),
     .          DISI2(NDOFC,*),       DISP2(NDOFC,*),
     .          VTANL(NDIME,NDIME-1,*)
C
      TWOPI=6.283185307179586D0
C
      NNODX=NNODL
      IF(NOCOL.EQ.1) THEN                          ! non-coincident mesh
       IF(ITASK.NE.1) NNODX=NNODN
      ENDIF
C
      IF(NMEMO1M.EQ.0) THEN
       DO IDIME=1,NDIME
        DO INODL=1,NNODL
         ELCO1(IDIME,INODL)=ELCOD(IDIME,INODL)
        ENDDO
       ENDDO
      ENDIF
C
      IF(NMEMO5M.EQ.0) THEN
       DO IDOFC=1,NDOFC
        DO INODL=1,NNODX
         ELDI1(IDOFC,INODL)=ELDIS(IDOFC,INODL)
        ENDDO
       ENDDO
      ENDIF
C
      IF(NMEMO2M.EQ.0) THEN
       IF(ITASK.EQ.1) THEN
        CALL SETM32(DVOLU,ELCO1,GPCOD,LNODS,PROPS,SHAPE,THICK,
     .              DERIV,POSGP,WEIGP,XJACM,VNORI,STRSG,VTANL)
       ELSE
        DO IGAUS=1,NGAUL
         DVOL1(IGAUS)=DVOLU(IGAUS)
         DO INODL=1,NNODX
          SHAP1(INODL,IGAUS)=SHAPE(INODL,IGAUS)
         ENDDO
         DO IDIME=1,NDIME
          GPCO1(IDIME,IGAUS)=GPCOD(IDIME,IGAUS)
          VNORL(IDIME,IGAUS)=VNORI(IDIME,IGAUS)
         ENDDO
        ENDDO
       ENDIF               ! itask.eq.1
      ELSE
C
C**** NMEMO2M=1 >> MECHANICAL PROBLEM INDEPENDENT OF ITASK
C
C     ITASK=1: checks negative volumes
C     ITASK>1: compute all the necessary arrays
C
       CALL SETM32(DVOL1,ELCO1,GPCO1,LNODS,PROPS,SHAP1,THICK,
     .             DERIV,POSGP,WEIGP,XJACM,VNORL,STRSG,VTANL)
      ENDIF                ! nmemo2m.eq.0
C
C**** UPDATES GEOMETRY FOR LARGE DISPLACEMENTS WITH NON-COINCID. MESH
C
      IF(NOCOL.EQ.1) THEN
       IF(LARGC.NE.0) THEN
        DO IDIME=1,NDIME
         DO INODE=1,NNODX
          ELDI1(IDIME,INODE)=ELCO1(IDIME,INODE)+ELDI1(IDIME,INODE)
         ENDDO
         DO INODE=1,NNODL
          ELCO2(IDIME,INODE)=ELDI1(IDIME,INODE)
          IF(LICOI.EQ.0) THEN    ! linearized computation of n
           IF(ITASK.EQ.4)        ! normal vector n is the same for R & K
     .      ELCO2(IDIME,INODE)=ELDI1(IDIME,INODE)-DISI2(IDIME,INODE)
           IF(ITASK.EQ.3.OR.ITASK.EQ.4)
     .      ELCO2(IDIME,INODE)=ELCO2(IDIME,INODE)-DISP2(IDIME,INODE)
          ENDIF
         END DO
        END DO
C
        CALL SETM32(DVOL2,ELCO2,GPCO2,LNODS,PROPS,SHAP2,THICK,
     .              DERIV,POSGP,WEIGP,XJACM,VNORL,STRSG,VTANL)
C
        IF(ITASK.EQ.1) THEN
         DO IGAUS=1,NGAUL
          DO IDIME=1,NDIME
           GPCO1(IDIME,IGAUS)=GPCO2(IDIME,IGAUS)
          ENDDO
         ENDDO
         DO IDIME=1,NDIME
          DO INODE=1,NNODL
           ELCO1(IDIME,INODE)=ELCO2(IDIME,INODE)
          END DO
         END DO
        ENDIF              ! itask.eq.1
        IF(ITASK.EQ.3.OR.ITASK.EQ.4) THEN
         DO IGAUS=1,NGAUL
          DVOL1(IGAUS)=DVOL2(IGAUS)
         ENDDO
        ENDIF              ! itask.eq.3.or.itask.eq.4
C
       ENDIF               ! largc.ne.0
      ENDIF                ! nocoi.eq.1
C
      RETURN
      END
