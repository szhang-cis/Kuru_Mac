      SUBROUTINE SMDG30(DVOLU,EMASS,SHAPE,STRSG,EHIST,
     .                  STREB,STREA,CTREB,CTREA,
     .                  CARTD,ELDIS,STRAN,XJACM,GPCOD,BSBAR)
C***********************************************************************
C
C**** THIS ROUTINE STORES NODAL STRESS & CONSTITUTIVE TENSOR FOR
C     DISCONTINUOUS GALERKIN METHOD
C
C     Note: NPRE4=NSTR1 & NPRE5=NKOST
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      COMMON/SMOGAU/FACTA
      COMMON/TUNING30/PPART,PPARB,PPARI,ESTAB
C
      DIMENSION DVOLU(*),       SHAPE(NNODS,*),
     .          STRSG(NSTR1,*), EMASS(NNODS,*)
      DIMENSION STREB(NSTR1,*), STREA(NNODS,*)
      DIMENSION CTREB(NPRE5,*), CTREA(NNODS,*), EHIST(NHIST,*)
      DIMENSION SAUXX(NSTR1),   CAUXX(NKOST)
      DIMENSION CARTD(NDIME,NNODS,*), ELDIS(NDIME,*)
      DIMENSION XJACM(NDIME,*),       XJACI(3,3)
      DIMENSION STRAN(NSTR1,*),
     .          SIGMA(3,3),           AUXIL(3,3)
      DIMENSION GPCOD(NDIME,*)
      DIMENSION BSBAR(NSTR1,NEVAB,*)
      DIMENSION XJACN(3,3),           XJANI(3,3)
      DIMENSION DMAT1(6,6),DMATT(3,3,3,3), AUXI1(3,3,3,3)
C
C**** COMPUTE TOTAL VOLUME IF NECESSARY
C
      IF(IABS(KSGAU).EQ.2) THEN
       FACTA=0.D+00
       DO IGAUS=1,NGAUL
        FACTA=FACTA+DVOLU(IGAUS)
       ENDDO
      ENDIF
C
C**** FOR 'ONE-POINT' INTEGRATION RULES
C
      IF(NGAUL.EQ.1) THEN
C$DIR SCALAR
       DO ISTRE=1,NSTR1
        DO INODE=1,NNODS
         STREB(ISTRE,INODE)=STRSG(ISTRE,1)
        ENDDO
       ENDDO
       DO ISTRE=1,NKOST
        DO INODE=1,NNODS
         CTREB(ISTRE,INODE)=EHIST(IPLAS(2)+ISTRE-1,1)
        ENDDO
       ENDDO
       RETURN
      ENDIF
C
C**** INITIALIZATION (not performed in elm000.f)
C
      DO ISTRE=1,NSTR1
       DO INODE=1,NNODS
        STREA(INODE,ISTRE)=0.0D0
        STREB(ISTRE,INODE)=0.0D0
       ENDDO
      ENDDO
      DO ISTRE=1,NKOST
       DO INODE=1,NNODS
        CTREA(INODE,ISTRE)=0.0D0
        CTREB(ISTRE,INODE)=0.0D0
       ENDDO
      ENDDO
C
C**** INTEGRAL PERFORMING
C
      KPARB=INT(PPARB)
C
C**** LOOP ON INTEGRATION POINTS
C
C$DIR SCALAR
      DO IGAUS=1,NGAUL
C
C**** CALCULATE THE JACOBIAN MATRIX
C
       DO ISTRE=1,NSTR1
        SAUXX(ISTRE)=STRSG(ISTRE,IGAUS)
       ENDDO
       DO ISTRE=1,NKOST
        CAUXX(ISTRE)=EHIST(IPLAS(2)+ISTRE-1,IGAUS)
       ENDDO
       IF(LARGE.EQ.1) THEN
        IF(KPARB.LE.0) THEN
         CALL LINEAR(CARTD(1,1,IGAUS), DEPSV,DESIG,DMATX,DSTRA,ELDIS,
     .               GPCOD(1,IGAUS),   LARGE,NDIME,NDOFN,NNODS,NSTR1,
     .               PROPS,            SHAPE(1,IGAUS),DUMMY,
     .               DUMMY,DUMMY,
     .               XJACM,XJACI,XJA3M,XJA3I,DETJM,DETJB,
     .               ELCOD,
     .                   1)
        ELSE
         CALL LINEAB(BSBAR(1,1,IGAUS), DEPSV,DESIG,DMATX,DSTRA,ELDIS,
     .               GPCOD(1,IGAUS),   LARGE,NDIME,NDOFN,NNODS,NSTR1,
     .               PROPS,            SHAPE(1,IGAUS),DUMMY,
     .               DUMMY,DUMMY,
     .               XJACM,XJACI,XJA3M,XJA3I,DETJM,DEFRA,
     .               ELCOD,CARTD(1,1,IGAUS),NEVAB,
     .               XJACN,XJANI,XJA3N,XJ3NI,DETJN,DETJB,
     .                   1)
        ENDIF
C
C**** BUILD STRESS TENSOR
C
        SIGMA(1,1)=SAUXX(1)
        IF(NDIME.EQ.3)THEN
         SIGMA(1,2)=SAUXX(3)
         SIGMA(1,3)=SAUXX(5)
         SIGMA(2,2)=SAUXX(2)
         SIGMA(2,3)=SAUXX(6)
         SIGMA(3,3)=SAUXX(4)
        ELSE IF(NDIME.EQ.2) THEN
         SIGMA(1,2)=SAUXX(3)
         SIGMA(2,2)=SAUXX(2)
        ENDIF
        DO 60 I=1,NDIME
        DO 60 J=I,NDIME
   60   SIGMA(J,I)=SIGMA(I,J)
C                                           
C**** COMPUTES FIRST PIOLA-KIRCHHOFF STRESS: P = F:S
C
        DO 70 I=1,NDIME
        DO 70 J=1,NDIME
        AUXIL(I,J)=0.D0
        DO 70 K=1,NDIME
   70   AUXIL(I,J)=AUXIL(I,J)+XJACM(I,K)*SIGMA(K,J)
        DO 80 I=1,NDIME
        DO 80 J=1,NDIME
   80   SIGMA(I,J)=AUXIL(I,J)  
C
        IF(NTYPE.EQ.1) THEN
         FPS33=DSQRT(2.0D0*STRAN(4,IGAUS)+1.0D0)
         DETJB=DETJB*FPS33
        ENDIF
C
        SAUXX(1)=SIGMA(1,1)/DETJB
        IF(NDIME.EQ.3)THEN
         SAUXX(3)=SIGMA(1,2)/DETJB
         SAUXX(5)=SIGMA(1,3)/DETJB
         SAUXX(2)=SIGMA(2,2)/DETJB
         SAUXX(6)=SIGMA(2,3)/DETJB
         SAUXX(4)=SIGMA(3,3)/DETJB
        ELSE IF(NDIME.EQ.2) THEN
         SAUXX(3)=SIGMA(1,2)/DETJB
         SAUXX(2)=SIGMA(2,2)/DETJB
        ENDIF
C                                                             T
C**** COMPUTES CONSTITUTIVE TENSOR DERIVED FROM P: C_P = F:C:F
C
        IKONT=0
        DO ISTRE=1,NSTRS
         INDEX=ISTRE
         IF(KSYMM.EQ.0) INDEX=1
         DO JSTRE=INDEX,NSTRS
          IKONT=IKONT+1
          DMAT1(ISTRE,JSTRE)=CAUXX(IKONT)
         ENDDO
        ENDDO
        IF(KSYMM.EQ.1) THEN
         DO ISTRE=1,NSTRS
          DO JSTRE=ISTRE,NSTRS
           DMAT1(JSTRE,ISTRE)=DMAT1(ISTRE,JSTRE)
          ENDDO
         ENDDO
        ENDIF
C
        DMATT(1,1,1,1)=DMAT1(1,1)
        IF(NDIME.NE.1) THEN
         DMATT(1,1,1,2)=DMAT1(1,3)
         DMATT(1,1,2,1)=DMAT1(1,3)
         DMATT(1,1,2,2)=DMAT1(1,2)
         DMATT(1,2,1,1)=DMAT1(3,1)
         DMATT(1,2,1,2)=DMAT1(3,3)
         DMATT(1,2,2,1)=DMAT1(3,3)
         DMATT(1,2,2,2)=DMAT1(3,2)
         DMATT(2,1,1,1)=DMAT1(3,1)
         DMATT(2,1,1,2)=DMAT1(3,3)
         DMATT(2,1,2,1)=DMAT1(3,3)
         DMATT(2,1,2,2)=DMAT1(3,2)
         DMATT(2,2,1,1)=DMAT1(2,1)
         DMATT(2,2,1,2)=DMAT1(2,3)
         DMATT(2,2,2,1)=DMAT1(2,3)
         DMATT(2,2,2,2)=DMAT1(2,2)
         IF(NDIME.EQ.3) THEN
          DMATT(1,1,3,3)=DMAT1(1,4)
          DMATT(1,2,3,3)=DMAT1(3,4)
          DMATT(2,1,3,3)=DMAT1(3,4)
          DMATT(2,2,3,3)=DMAT1(2,4)
          DMATT(3,3,1,1)=DMAT1(4,1)
          DMATT(3,3,1,2)=DMAT1(4,3)
          DMATT(3,3,2,1)=DMAT1(4,3)
          DMATT(3,3,2,2)=DMAT1(4,2)
          DMATT(3,3,3,3)=DMAT1(4,4)
          DMATT(1,1,1,3)=DMAT1(1,5)
          DMATT(1,1,2,3)=DMAT1(1,6)
          DMATT(1,1,3,1)=DMAT1(1,5)
          DMATT(1,1,3,2)=DMAT1(1,6)
          DMATT(1,2,1,3)=DMAT1(3,5)
          DMATT(1,2,2,3)=DMAT1(3,6)
          DMATT(1,2,3,1)=DMAT1(3,5)
          DMATT(1,2,3,2)=DMAT1(3,6)
          DMATT(1,3,1,1)=DMAT1(5,1)
          DMATT(1,3,1,2)=DMAT1(5,3)
          DMATT(1,3,1,3)=DMAT1(5,5)
          DMATT(1,3,2,1)=DMAT1(5,3)
          DMATT(1,3,2,2)=DMAT1(5,2)
          DMATT(1,3,2,3)=DMAT1(5,6)
          DMATT(1,3,3,1)=DMAT1(5,5)
          DMATT(1,3,3,2)=DMAT1(5,6)
          DMATT(1,3,3,3)=DMAT1(5,4)
          DMATT(2,1,1,3)=DMAT1(3,5)
          DMATT(2,1,2,3)=DMAT1(3,6)
          DMATT(2,1,3,1)=DMAT1(3,5)
          DMATT(2,1,3,2)=DMAT1(3,6)
          DMATT(2,2,1,3)=DMAT1(2,5)
          DMATT(2,2,2,3)=DMAT1(2,6)
          DMATT(2,2,3,1)=DMAT1(2,5)
          DMATT(2,2,3,2)=DMAT1(2,6)
          DMATT(2,3,1,1)=DMAT1(6,1)
          DMATT(2,3,1,2)=DMAT1(6,3)
          DMATT(2,3,1,3)=DMAT1(6,5)
          DMATT(2,3,2,1)=DMAT1(6,3)
          DMATT(2,3,2,2)=DMAT1(6,2)
          DMATT(2,3,2,3)=DMAT1(6,6)
          DMATT(2,3,3,1)=DMAT1(6,5)
          DMATT(2,3,3,2)=DMAT1(6,6)
          DMATT(2,3,3,3)=DMAT1(6,4)
          DMATT(3,1,1,1)=DMAT1(5,1)
          DMATT(3,1,1,2)=DMAT1(5,3)
          DMATT(3,1,1,3)=DMAT1(5,5)
          DMATT(3,1,2,1)=DMAT1(5,3)
          DMATT(3,1,2,2)=DMAT1(5,2)
          DMATT(3,1,2,3)=DMAT1(5,6)
          DMATT(3,1,3,1)=DMAT1(5,5)
          DMATT(3,1,3,2)=DMAT1(5,6)
          DMATT(3,1,3,3)=DMAT1(5,4)
          DMATT(3,2,1,1)=DMAT1(6,1)
          DMATT(3,2,1,2)=DMAT1(6,3)
          DMATT(3,2,1,3)=DMAT1(6,5)
          DMATT(3,2,2,1)=DMAT1(6,3)
          DMATT(3,2,2,2)=DMAT1(6,2)
          DMATT(3,2,2,3)=DMAT1(6,6)
          DMATT(3,2,3,1)=DMAT1(6,5)
          DMATT(3,2,3,2)=DMAT1(6,6)
          DMATT(3,2,3,3)=DMAT1(6,4)
          DMATT(3,3,1,3)=DMAT1(4,5)
          DMATT(3,3,2,3)=DMAT1(4,6)
          DMATT(3,3,3,1)=DMAT1(4,5)
          DMATT(3,3,3,2)=DMAT1(4,6)
         ENDIF                     ! ndime.eq.3
        ENDIF                      ! ndime.ne.1
C
C**** COMPUTE CONSTITUTIVE TENSOR C_P
C
        DO 170 I=1,NDIME
        DO 170 J=1,NDIME
        DO 170 K=1,NDIME
        DO 170 L=1,NDIME
        AUXI1(I,J,K,L)=0.D0
        DO 170 I1=1,NDIME
        DO 170 L1=1,NDIME
  170   AUXI1(I,J,K,L)=AUXI1(I,J,K,L)+
     .                 XJACM(I,I1)*DMATT(I1,J,K,L1)*XJACM(L,L1)
C
        DO 180 I=1,NDIME
        DO 180 J=1,NDIME
        DO 180 K=1,NDIME
        DO 180 L=1,NDIME
  180   DMATT(I,J,K,L)=AUXI1(I,J,K,L)
C
        DMAT1(1,1)=DMATT(1,1,1,1)
        IF(NDIME.NE.1) THEN
         DMAT1(1,3)=DMATT(1,1,1,2)
         DMAT1(1,3)=DMATT(1,1,2,1)
         DMAT1(1,2)=DMATT(1,1,2,2)
         DMAT1(3,1)=DMATT(1,2,1,1)
         DMAT1(3,3)=DMATT(1,2,1,2)
         DMAT1(3,3)=DMATT(1,2,2,1)
         DMAT1(3,2)=DMATT(1,2,2,2)
         DMAT1(3,1)=DMATT(2,1,1,1)
         DMAT1(3,3)=DMATT(2,1,1,2)
         DMAT1(3,3)=DMATT(2,1,2,1)
         DMAT1(3,2)=DMATT(2,1,2,2)
         DMAT1(2,1)=DMATT(2,2,1,1)
         DMAT1(2,3)=DMATT(2,2,1,2)
         DMAT1(2,3)=DMATT(2,2,2,1)
         DMAT1(2,2)=DMATT(2,2,2,2)
         IF(NDIME.EQ.3) THEN
          DMAT1(1,4)=DMATT(1,1,3,3)
          DMAT1(3,4)=DMATT(1,2,3,3)
          DMAT1(3,4)=DMATT(2,1,3,3)
          DMAT1(2,4)=DMATT(2,2,3,3)
          DMAT1(4,1)=DMATT(3,3,1,1)
          DMAT1(4,3)=DMATT(3,3,1,2)
          DMAT1(4,3)=DMATT(3,3,2,1)
          DMAT1(4,2)=DMATT(3,3,2,2)
          DMAT1(4,4)=DMATT(3,3,3,3)
          DMAT1(1,5)=DMATT(1,1,1,3)
          DMAT1(1,6)=DMATT(1,1,2,3)
          DMAT1(3,5)=DMATT(1,2,1,3)
          DMAT1(3,6)=DMATT(1,2,2,3)
          DMAT1(1,5)=DMATT(1,1,3,1)
          DMAT1(1,6)=DMATT(1,1,3,2)
          DMAT1(3,5)=DMATT(1,2,3,1)
          DMAT1(3,6)=DMATT(1,2,3,2)
          DMAT1(5,1)=DMATT(1,3,1,1)
          DMAT1(5,3)=DMATT(1,3,1,2)
          DMAT1(5,5)=DMATT(1,3,1,3)
          DMAT1(5,3)=DMATT(1,3,2,1)
          DMAT1(5,2)=DMATT(1,3,2,2)
          DMAT1(5,6)=DMATT(1,3,2,3)
          DMAT1(5,5)=DMATT(1,3,3,1)
          DMAT1(5,6)=DMATT(1,3,3,2)
          DMAT1(5,4)=DMATT(1,3,3,3)
          DMAT1(3,5)=DMATT(2,1,1,3)
          DMAT1(3,6)=DMATT(2,1,2,3)
          DMAT1(3,5)=DMATT(2,1,3,1)
          DMAT1(3,6)=DMATT(2,1,3,2)
          DMAT1(2,5)=DMATT(2,2,1,3)
          DMAT1(2,6)=DMATT(2,2,2,3)
          DMAT1(2,5)=DMATT(2,2,3,1)
          DMAT1(2,6)=DMATT(2,2,3,2)
          DMAT1(6,1)=DMATT(2,3,1,1)
          DMAT1(6,3)=DMATT(2,3,1,2)
          DMAT1(6,5)=DMATT(2,3,1,3)
          DMAT1(6,3)=DMATT(2,3,2,1)
          DMAT1(6,2)=DMATT(2,3,2,2)
          DMAT1(6,6)=DMATT(2,3,2,3)
          DMAT1(6,5)=DMATT(2,3,3,1)
          DMAT1(6,6)=DMATT(2,3,3,2)
          DMAT1(6,4)=DMATT(2,3,3,3)
          DMAT1(5,1)=DMATT(3,1,1,1)
          DMAT1(5,3)=DMATT(3,1,1,2)
          DMAT1(5,5)=DMATT(3,1,1,3)
          DMAT1(5,3)=DMATT(3,1,2,1)
          DMAT1(5,2)=DMATT(3,1,2,2)
          DMAT1(5,6)=DMATT(3,1,2,3)
          DMAT1(5,5)=DMATT(3,1,3,1)
          DMAT1(5,6)=DMATT(3,1,3,2)
          DMAT1(5,4)=DMATT(3,1,3,3)
          DMAT1(6,1)=DMATT(3,2,1,1)
          DMAT1(6,3)=DMATT(3,2,1,2)
          DMAT1(6,5)=DMATT(3,2,1,3)
          DMAT1(6,3)=DMATT(3,2,2,1)
          DMAT1(6,2)=DMATT(3,2,2,2)
          DMAT1(6,6)=DMATT(3,2,2,3)
          DMAT1(6,5)=DMATT(3,2,3,1)
          DMAT1(6,6)=DMATT(3,2,3,2)
          DMAT1(6,4)=DMATT(3,2,3,3)
          DMAT1(4,5)=DMATT(3,3,1,3)
          DMAT1(4,6)=DMATT(3,3,2,3)
          DMAT1(4,5)=DMATT(3,3,3,1)
          DMAT1(4,6)=DMATT(3,3,3,2)
         ENDIF                     ! ndime.eq.3
        ENDIF                      ! ndime.ne.1
C
        IKONT=0
        DO ISTRS=1,NSTRS
         INDEX=ISTRS
         IF(KSYMM.EQ.0) INDEX=1
         DO JSTRS=INDEX,NSTRS
          IKONT=IKONT+1
          CAUXX(IKONT)=DMAT1(ISTRS,JSTRS)
         ENDDO
        ENDDO
       ENDIF             ! large.eq.1
C
       DO ISTRE=1,NSTR1
        DO INODE=1,NNODS
         STREA(INODE,ISTRE)=STREA(INODE,ISTRE)+SHAPE(INODE,IGAUS)*
     .                                         SAUXX(ISTRE)*DVOLU(IGAUS)
        ENDDO
       ENDDO
       DO ISTRE=1,NKOST
        DO INODE=1,NNODS
         CTREA(INODE,ISTRE)=CTREA(INODE,ISTRE)+SHAPE(INODE,IGAUS)*
     .                                         CAUXX(ISTRE)*DVOLU(IGAUS)
        ENDDO
       ENDDO
      ENDDO              ! igaus=1,ngaul
C
C**** PERFORMS MATRIX PRODUCT
C
C$DIR SCALAR
      DO ISTRE=1,NSTR1
       DO INODE=1,NNODS
        DO JNODE=1,NNODS
         STREB(ISTRE,INODE)=STREB(ISTRE,INODE)+
     .                      STREA(JNODE,ISTRE)*EMASS(INODE,JNODE)
        ENDDO
       ENDDO
      ENDDO
      DO ISTRE=1,NKOST
       DO INODE=1,NNODS
        DO JNODE=1,NNODS
         CTREB(ISTRE,INODE)=CTREB(ISTRE,INODE)+
     .                      CTREA(JNODE,ISTRE)*EMASS(INODE,JNODE)
        ENDDO
       ENDDO
      ENDDO
C
      RETURN
      END
