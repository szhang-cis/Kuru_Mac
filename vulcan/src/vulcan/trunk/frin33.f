      SUBROUTINE FRIN33(PROPS,LNODS,ELDIS,BMSIG,DVOLU,SHAPE,EHIST,TENOD,
     .                  ELCOD,PWOEL,PREAL,TGAPL,VNORL,VTANL,DISIL,DISPL,
     .                  STDGL,CTDGL)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE CONTACT FORCES FOR THE DG ELEMENT
C     ( ELEMENT NO. 33 )
C
C     Notes:
C     For DG linking problems, NDIME=NDOFN
C
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
      DIMENSION PROPS(*),       ELDIS(NDOFC,*),
     .          BMSIG(*),       DVOLU(*),
     .          SHAPE(NNODN,*), EHIST(NHIST,*)
      DIMENSION TENOD(*),       ELCOD(NDIME,*),
     .          PWOEL(*),       PREAL(*),
     .          TGAPL(*),       LNODS(*)
      DIMENSION VNORL(NDIME,*), DGAPS(3),
     .          VTANL(NDIME,NDIME-1,*),
     .          DISIL(NDOFC,*), DISPL(NDOFC,*)
      DIMENSION DGAPT(2),       PREST(2),
     .          RIGI1(3,3),     RIGI2(2,3), RIGI3(2,2)
      DIMENSION STDGL(NPRE4,*), CTDGL(NPRE5,*), TDG(3)
      DIMENSION STRE1(6),   STRE2(6),  SIGM1(3,3), SIGM2(3,3)
      DIMENSION CTRE1(36),  CTRE2(36), DMAT1(6,6), DMAT2(6,6),
     .          VNOR1(6),   VAUXN(6),  DMATX(3,3)
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
C**** INITIALISE MECHANICAL COUPLING TERM
C
      IF(ITERME.GT.0) THEN              ! bidirectional coupling
       DO INODL=1,NNODL
        PWOEL(INODL)=0.0D0
       END DO
      ENDIF
C
C**** LOOP OVER INTEGRATION POINTS
C
      DO 101 IGAUS=1,NGAUL
C
C**** COMPUTES TEMPERATURE
C
      TGAUS=0.0D0
      IF(ITERME.GE.0) THEN              ! coupled problems
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
C**** COMPUTES THE FIRST FORCE VECTOR-TYPE TERM (T_DG)
C
C     Note: for simplicity, the stiffness matrix corresponding to this
C           term is neglected (see mass33.f)
C
      DO ISTRE=1,NSTR1
       STRE1(ISTRE)=0.0D0               ! stress at boundary 1
       STRE2(ISTRE)=0.0D0               ! stress at boundary 2
       DO INODL=1,NNOBO
        STRE1(ISTRE)=STRE1(ISTRE)+SHAPE(INODL,IGAUS)*STDGL(ISTRE,INODL)
        STRE2(ISTRE)=STRE2(ISTRE)+
     .                       SHAPE(INODL,IGAUS)*STDGL(ISTRE,INODL+NNOBO)
       END DO
      END DO
C
      SIGM1(1,1)=STRE1(1)               ! build stress tensor at b1
      IF(NDIME.EQ.3) THEN
       SIGM1(1,2)=STRE1(3)
       SIGM1(1,3)=STRE1(5)
       SIGM1(2,2)=STRE1(2)
       SIGM1(2,3)=STRE1(6)
       SIGM1(3,3)=STRE1(4)
      ELSE IF(NDIME.EQ.2) THEN
       SIGM1(1,2)=STRE1(3)
       SIGM1(2,2)=STRE1(2)
      ENDIF
C
      SIGM2(1,1)=STRE2(1)               ! build stress tensor at b2
      IF(NDIME.EQ.3) THEN
       SIGM2(1,2)=STRE2(3)
       SIGM2(1,3)=STRE2(5)
       SIGM2(2,2)=STRE2(2)
       SIGM2(2,3)=STRE2(6)
       SIGM2(3,3)=STRE2(4)
      ELSE IF(NDIME.EQ.2) THEN
       SIGM2(1,2)=STRE2(3)
       SIGM2(2,2)=STRE2(2)
      ENDIF
C
      DO 60 I=1,NDIME
      DO 60 J=I,NDIME
      SIGM1(J,I)=SIGM1(I,J)
   60 SIGM2(J,I)=SIGM2(I,J)
C
      DO IDIME=1,NDIME                  ! computes <P>.N product
       TDG(IDIME)=0.0D0
       DO JDIME=1,NDIME
        TDG(IDIME)=TDG(IDIME)+
     .  (SIGM1(IDIME,JDIME)+SIGM2(IDIME,JDIME))/2.0D0*VNORL(JDIME,IGAUS)
       ENDDO
      ENDDO
C
      DO INODE=1,NNOBO
       SAREA=SHAPE(INODE,IGAUS)*DVOLU(IGAUS)
       IEVAB=(INODE-1)*NDIME
C$DIR SCALAR
       DO IDIME=1,NDIME
        IEVAB=IEVAB+1
        BMSIG(IEVAB)=BMSIG(IEVAB)+TDG(IDIME)*SAREA
       ENDDO
      ENDDO
C
C**** COMPUTES THE SECOND MASS MATRIX-TYPE TERM (C_DG)
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
      DO IDOFN=1,NDOFN                  ! computes gap (NDOFN=NDIME)
       DGAPS(IDOFN)=0.0D+00
       DO INODL=1,NNOBO
        ELDI1=ELDIS(IDOFN,INODL)
        ELDI2=ELDIS(IDOFN,INODL+NNOBO)
        DGAPS(IDOFN)=DGAPS(IDOFN)+SHAPE(INODL,IGAUS)*(ELDI1-ELDI2)
       ENDDO
      ENDDO
C
      DO IDIME=1,NDIME                  ! computes (C:NxN)*gDG product
       TDG(IDIME)=0.0D0
       DO JDIME=1,NDIME
        TDG(IDIME)=TDG(IDIME)+DMATX(IDIME,JDIME)*DGAPS(JDIME)
       ENDDO
      ENDDO
C
      DO INODE=1,NNOBO
       SAREA=SHAPE(INODE,IGAUS)*DVOLU(IGAUS)
       IEVAB=(INODE-1)*NDIME
C$DIR SCALAR
       DO IDIME=1,NDIME
        IEVAB=IEVAB+1
        BMSIG(IEVAB)=BMSIG(IEVAB)+TDG(IDIME)*SAREA
       ENDDO
      ENDDO
C
  101 CONTINUE
C
      DO IEVAB=1,NEVBO
       BMSIG(IEVAB+NEVBO)=-BMSIG(IEVAB)
      END DO
C
      RETURN
      END
