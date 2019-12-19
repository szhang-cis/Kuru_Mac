      SUBROUTINE PUFOBA(STRSG,STRSGL,XJACM,XJA3M,DETJM,IFOBA)
C***********************************************************************
C
C**** THIS ROUTINE PERFORMS PUSH-FORWARD & PULL-BACK OPERATIONS ON:
C
C     IFOBA=1 > STRESS TENSOR
C     IFOBA=2 > STRAIN TENSOR
C     IFOBA=3 > NO TRANSFORMATION
C
C     Notes: STRSG(4) (strains, stresses, etc.) must be transformed for
C            plane stress problems (NTYPE=1)
C
C     Input: STRSG, Output: STRSGL
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION STRSG(*),       STRSGL(*)
      DIMENSION XJACM(NDIME,*)
      DIMENSION SIGMA(3,3),     AUXIL(3,3)
C
      IF(IFOBA.EQ.3) THEN                                    ! no change
       DO ISTR1=1,NSTR1
        STRSGL(ISTR1)=STRSG(ISTR1)
       ENDDO
       RETURN
      ENDIF
C
C**** PUSH-FORWARD/PULL-BACK OPERATIONS
C
      SIGMA(1,1)=STRSG(1)
      IF(NTYPE.EQ.4)THEN
       SIGMA(1,2)=STRSG(3)
       SIGMA(1,3)=STRSG(5)
       SIGMA(2,2)=STRSG(2)
       SIGMA(2,3)=STRSG(6)
       SIGMA(3,3)=STRSG(4)
       IF(IFOBA.EQ.2) THEN
        SIGMA(1,2)=SIGMA(1,2)/2.0D0
        SIGMA(1,3)=SIGMA(1,3)/2.0D0
        SIGMA(2,3)=SIGMA(2,3)/2.0D0
       ENDIF
      ELSE IF(NTYPE.NE.5) THEN
       SIGMA(1,2)=STRSG(3)
       SIGMA(2,2)=STRSG(2)
       SIGMA(3,3)=STRSG(4)
       IF(IFOBA.EQ.2) SIGMA(1,2)=SIGMA(1,2)/2.0D0
      ENDIF
      DO 60 I=1,NDIME
      DO 60 J=I,NDIME
   60 SIGMA(J,I)=SIGMA(I,J)
C
      DO 70 I=1,NDIME
      DO 70 J=1,NDIME
      AUXIL(I,J)=0.0D0
      DO 70 K=1,NDIME
      DO 70 L=1,NDIME
      IF(IFOBA.EQ.1) THEN
       AUXIL(I,J)=AUXIL(I,J)+XJACM(I,K)*SIGMA(K,L)*XJACM(J,L)
      ELSE           ! ifoba=2
       AUXIL(I,J)=AUXIL(I,J)+XJACM(K,I)*SIGMA(K,L)*XJACM(L,J)
      ENDIF
   70 CONTINUE
      DO 80 I=1,NDIME
      DO 80 J=1,NDIME
   80 SIGMA(I,J)=AUXIL(I,J)  
      IF(NTYPE.EQ.3) SIGMA(3,3)=XJA3M*SIGMA(3,3)*XJA3M
C
      STRSGL(1)=SIGMA(1,1)/DETJM
      IF(NTYPE.EQ.4)THEN
       STRSGL(3)=SIGMA(1,2)/DETJM
       STRSGL(5)=SIGMA(1,3)/DETJM
       STRSGL(2)=SIGMA(2,2)/DETJM
       STRSGL(6)=SIGMA(2,3)/DETJM
       STRSGL(4)=SIGMA(3,3)/DETJM
       IF(IFOBA.EQ.2) THEN
        STRSGL(3)=STRSGL(3)*2.0D0
        STRSGL(5)=STRSGL(5)*2.0D0
        STRSGL(6)=STRSGL(6)*2.0D0
       ENDIF
      ELSE IF(NTYPE.NE.5) THEN
       STRSGL(3)=SIGMA(1,2)/DETJM
       STRSGL(2)=SIGMA(2,2)/DETJM
       STRSGL(4)=SIGMA(3,3)/DETJM
       IF(IFOBA.EQ.2) STRSGL(3)=STRSGL(3)*2.0D0
      ENDIF
C
      RETURN
      END
