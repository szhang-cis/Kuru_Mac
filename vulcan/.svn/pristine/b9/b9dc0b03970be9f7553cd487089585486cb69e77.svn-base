      SUBROUTINE MASS33(DVOLU,PROPS,LNODS,SHAPE,WSTIF,EHIST,ELDIS,TENOD,
     .                  DISPL,ELCOD,VNORL,DMATX,AUXS1,AUXS2,BMATX,VTANL,
     .                  CTDGL)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE STIFFNESS MATRIX FOR THE DG ELEMENT
C     ( ELEMENT NO. 33 )
C
C     Note:
C     SHAPE must be dimensioned with NNOBO
C     ELDIS & TENOD must be dimensioned with NNODL
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
      COMMON/TUNING4/RITEN,RITEF,TRUPL,TRUPM,TOLGA,TOLGAM
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
      DIMENSION CTDGL(NPRE5,*)
      DIMENSION CTRE1(36),  CTRE2(36), DMAT1(6,6), DMAT2(6,6),
     .          VNOR1(6),   VAUXN(6)
C
      NNOBO=NNODL/2
      NEVBO=NNOBO*NDOFN
C
C**** PROPS(2)=BETA DG PARAMETER
C     PROPS(5)=LIQUIDUS TEMPERATURE
C
      BETDG=PROPS(2)
      TEMPL=PROPS(5)
      INOTE=INT(PROPS(6))
C
C**** CHARACTERISTIC LENGTH
C
      HACDG=0.0D0
      DO IDIME=1,NDIME
       HACDG=HACDG+(ELCOD(IDIME,1)*ELCOD(IDIME,1)-
     .              ELCOD(IDIME,2)*ELCOD(IDIME,2))
      ENDDO
      HACDG=DSQRT(HACDG)
C
C**** LOOP ON INTEGRATION POINTS
C
      DO 101 IGAUS=1,NGAUL
C
C**** COMPUTES TEMPERATURE
C
      TGAUS=0.0D0
      IF(ITERME.GE.0) THEN                 ! coupled problems
       IF(INOTE.EQ.1) THEN
        DO INODL=1,NNOBO
         TGAUS=TGAUS+SHAPE(INODL,IGAUS)*TENOD(INODL)
        END DO
       ELSE
        DO INODL=1,NNOBO
         TGAUS=TGAUS+SHAPE(INODL,IGAUS)*TENOD(INODL+NNOBO)
        END DO
       ENDIF
      ENDIF
C
C**** COMPUTES THE SECOND MASS MATRIX-TYPE TERM (C_DG)
C
C     Note: the first mass/stiffness matrix-type term (A_DG) is not
C           computed (see frin33.f)
C
      DO ISTRE=1,NKOST
       CTRE1(ISTRE)=0.0D0               ! constitutive tensor at b1
       CTRE2(ISTRE)=0.0D0               ! constitutive tensor at b2
       DO INODL=1,NNOBO
        CTRE1(ISTRE)=CTRE1(ISTRE)+SHAPE(INODL,IGAUS)*CTDGL(ISTRE,INODL)
        CTRE2(ISTRE)=CTRE2(ISTRE)+
     .                       SHAPE(INODL,IGAUS)*CTDGL(ISTRE,INODL+NNOBO)
       END DO
      END DO
C
      IKONT=0
      DO ISTRE=1,NSTRS
       INDEX=ISTRE
       IF(KSYMM.EQ.0) INDEX=1
       DO JSTRE=INDEX,NSTRS
        IKONT=IKONT+1
        DMAT1(ISTRE,JSTRE)=CTRE1(IKONT) ! build stress tensor at b1
        DMAT2(ISTRE,JSTRE)=CTRE2(IKONT) ! build stress tensor at b1
       ENDDO
      ENDDO
      IF(KSYMM.EQ.1) THEN
       DO ISTRE=1,NSTRS
        DO JSTRE=ISTRE,NSTRS
         DMAT1(JSTRE,ISTRE)=DMAT1(ISTRE,JSTRE)
         DMAT2(JSTRE,ISTRE)=DMAT2(ISTRE,JSTRE)
        ENDDO
       ENDDO
      ENDIF
C
      VNOR1(1)=VNORL(1,IGAUS)*VNORL(1,IGAUS) ! NxN product
      IF(NDIME.EQ.3) THEN
       VNOR1(2)=VNORL(2,IGAUS)*VNORL(2,IGAUS)
       VNOR1(3)=VNORL(1,IGAUS)*VNORL(2,IGAUS)*2.0D0
       VNOR1(4)=VNORL(3,IGAUS)*VNORL(3,IGAUS)
       VNOR1(5)=VNORL(1,IGAUS)*VNORL(3,IGAUS)*2.0D0
       VNOR1(6)=VNORL(2,IGAUS)*VNORL(3,IGAUS)*2.0D0
      ELSE IF(NDIME.EQ.2) THEN
       VNOR1(2)=VNORL(2,IGAUS)*VNORL(2,IGAUS)
       VNOR1(3)=VNORL(1,IGAUS)*VNORL(2,IGAUS)*2.0D0
      ENDIF
C
      DO ISTRE=1,NSTRS                  ! computes <C_P>:NxN product
       VAUXN(ISTRE)=0.0D0
       DO JSTRE=1,NSTRS
        VAUXN(ISTRE)=VAUXN(ISTRE)+BETDG/HACDG*
     .        (DMAT1(ISTRE,JSTRE)+DMAT2(ISTRE,JSTRE))/2.0D0*VNOR1(JSTRE)
       ENDDO
      ENDDO
      DMATX(1,1)=VAUXN(1)
      IF(NDIME.EQ.3) THEN
       DMATX(1,2)=VAUXN(3)
       DMATX(1,3)=VAUXN(5)
       DMATX(2,2)=VAUXN(2)
       DMATX(2,3)=VAUXN(6)
       DMATX(3,3)=VAUXN(4)
      ELSE IF(NDIME.EQ.2) THEN
       DMATX(1,2)=VAUXN(3)
       DMATX(2,2)=VAUXN(2)
      ENDIF
      DO IDIME=1,NDIME
       DO JDIME=IDIME,NDIME
        DMATX(JDIME,IDIME)=DMATX(IDIME,JDIME)
       ENDDO
      ENDDO
C
C**** COMPUTES CONTACT CONSTITUTIVE TENSOR
C
      CALL WBOUND(NDIME,NEVBO,NNOBO,
     .            SHAPE(1,IGAUS),WSTIF,DMATX,BMATX,
     .            KSYMM,AUXS1)
C
      CALL WINTERT(DVOLU(IGAUS),NEVAB,NEVBO,WSTIF,KSYMM,
     .             AUXS1,AUXS2)
C
  101 CONTINUE
C
      RETURN
      END
