      SUBROUTINE BBMEAN(CARTD,CMEAN,GPCOD,SHAPE,DVOLU,NAUXI)
C***********************************************************************
C
C**** THIS SUBROUTINE EVALUATES THE VOLUMETRIC-BAR STRAIN DISPLACEMENT
C     MATRIX CMEAN(3,NNODL)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION CARTD(NDIME,NNODL,*), CMEAN(NAUXI,*),
     .          GPCOD(NDIME,*),       SHAPE(NNODL,*),
     .          DVOLU(*)
C
      GO TO(1,2,3) NDIME
C
    1 CONTINUE
C
C**** 1D CASE (not implemented yet)
C
      CALL RUNEND('ERROR IN BBMEAN; NDIME = 1         ')
      RETURN
C
    2 CONTINUE
C
C**** 2D CASE FOR 3-NODED & 4-NODED ELEMENTS
C
      ICHEQ=0
      IF(NNODL.EQ.3.AND.NGAUL.EQ.1) THEN
       ICHEQ=1
       VOLIN=1.0/DVOLU(1)
      ENDIF
C
      IF(NNODL.EQ.4.AND.NGAUL.EQ.1) THEN
       ICHEQ=1
       VOLIN=1.0/DVOLU(1)
      ENDIF
C
      IF(NNODL.EQ.4.AND.NGAUL.EQ.3) THEN
       ICHEQ=1
       VOLIN=1.0/(DVOLU(1)+DVOLU(2)+DVOLU(3))
      ENDIF
C
      IF(NNODL.EQ.4.AND.NGAUL.EQ.6) THEN
       ICHEQ=1
       VOLIN=1.0/(DVOLU(1)+DVOLU(2)+DVOLU(3)+DVOLU(4)+DVOLU(5)+DVOLU(6))
      ENDIF
C
      IF(NNODL.EQ.4.AND.NGAUL.EQ.7) THEN
       ICHEQ=1
       VOLIN=1.0/(DVOLU(1)+DVOLU(2)+DVOLU(3)+DVOLU(4)+DVOLU(5)+DVOLU(6)+
     .            DVOLU(7))
      ENDIF
C
      IF(NNODL.EQ.4.AND.NGAUL.EQ.4) THEN
       ICHEQ=1
       VOLIN=1.0/(DVOLU(1)+DVOLU(2)+DVOLU(3)+DVOLU(4))
      ENDIF
C
C**** APPROXIMATED FORM FOR 8-NODED & 9-NODED ELEMENTS
C
      IF((NNODL.EQ.8.OR.NNODL.EQ.9).AND.NGAUL.EQ.4) THEN
       ICHEQ=1
       VOLIN=1.0/(DVOLU(1)+DVOLU(2)+DVOLU(3)+DVOLU(4))
      ENDIF
C
      IF((NNODL.EQ.8.OR.NNODL.EQ.9).AND.NGAUL.EQ.9) THEN
       ICHEQ=1
       VOLIN=1.0/(DVOLU(1)+DVOLU(2)+DVOLU(3)+DVOLU(4)+
     .            DVOLU(5)+DVOLU(6)+DVOLU(7)+DVOLU(8)+DVOLU(9))
      ENDIF
C
      IF(ICHEQ.EQ.0)
     . CALL RUNEND('ERROR IN BBMEAN; NDIME = 2         ')
C
      DO 20 INODE=1,NNODL
      CMEAN(1,INODE)=0.0D0
      CMEAN(2,INODE)=0.0D0
      CMEAN(3,INODE)=0.0D0
C
      DO 10 IGAUS=1,NGAUL
      CMEAN(1,INODE)=CMEAN(1,INODE)+VOLIN*CARTD(1,INODE,IGAUS)*
     .                              DVOLU(IGAUS)
      CMEAN(2,INODE)=CMEAN(2,INODE)+VOLIN*CARTD(2,INODE,IGAUS)*
     .                              DVOLU(IGAUS)
      IF(NTYPE.EQ.3) THEN
       GPCOA=0.0D0                                      ! average radius
       DO IGAUX=1,NGAUL
        GPCOA=GPCOA+GPCOD(1,IGAUX)
       ENDDO
       GPCOA=GPCOA/NGAUL
       IF(GPCOD(1,IGAUS).GT.(1.0D-10*GPCOA)) THEN     ! see GPCOA in *.f
        CMEAN(3,INODE)=CMEAN(3,INODE)+VOLIN*SHAPE(INODE,IGAUS)/
     .                                GPCOD(1,IGAUS)*DVOLU(IGAUS)
       ENDIF
      ENDIF
   10 CONTINUE     
   20 CONTINUE     
      RETURN
C
    3 CONTINUE
C
C**** 3D CASE FOR 4-NODED & 8-NODED ELEMENTS
C
      ICHEQ=0
      IF(NNODL.EQ.4.AND.NGAUL.EQ.1) THEN
       ICHEQ=1
       VOLIN=1.0/DVOLU(1)
      ENDIF
C
      IF(NNODL.EQ.5.AND.NGAUL.EQ.4) THEN
       ICHEQ=1
       VOLIN=1.0/(DVOLU(1)+DVOLU(2)+DVOLU(3)+DVOLU(4))
      ENDIF
C
      IF(NNODL.EQ.5.AND.NGAUL.EQ.5) THEN
       ICHEQ=1
       VOLIN=1.0/(DVOLU(1)+DVOLU(2)+DVOLU(3)+DVOLU(4)+DVOLU(5))
      ENDIF
C
      IF(NNODL.EQ.8.AND.NGAUL.EQ.8) THEN
       ICHEQ=1
       VOLIN=1.0/(DVOLU(1)+DVOLU(2)+DVOLU(3)+DVOLU(4)+
     .            DVOLU(5)+DVOLU(6)+DVOLU(7)+DVOLU(8))
      ENDIF
C
      IF(ICHEQ.EQ.0)
     . CALL RUNEND('ERROR IN BBMEAN; NDIME = 3         ')
C
      DO 40 INODE=1,NNODL
      CMEAN(1,INODE)=0.0D0
      CMEAN(2,INODE)=0.0D0
      CMEAN(3,INODE)=0.0D0
C
      DO 30 IGAUS=1,NGAUL
      CMEAN(1,INODE)=CMEAN(1,INODE)+VOLIN*CARTD(1,INODE,IGAUS)*
     .                              DVOLU(IGAUS)
      CMEAN(2,INODE)=CMEAN(2,INODE)+VOLIN*CARTD(2,INODE,IGAUS)*
     .                              DVOLU(IGAUS)
      CMEAN(3,INODE)=CMEAN(3,INODE)+VOLIN*CARTD(3,INODE,IGAUS)*
     .                              DVOLU(IGAUS)
   30 CONTINUE     
   40 CONTINUE     
      RETURN
      END
