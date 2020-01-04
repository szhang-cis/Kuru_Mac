      SUBROUTINE BBMATX(BSBAR,CMEAN,CARTD,GPCOD,SHAPE)
C***********************************************************************
C
C**** THIS SUBROUTINE EVALUATES THE B-BAR STRAIN-DISPLACEMENT MATRIX
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION BSBAR(NSTR1,*), CARTD(NDIME,*),
     .          CMEAN(3,*),     GPCOD(*),
     .          SHAPE(*)
C
c     DATA C1/0.33333333333/,C2/0.66666666666/
      C1=0.33333333333D0
      C2=0.66666666666D0
C
      NGASH=0
C
      DO INODE=1,NNODL
       TEMP1=C1*CMEAN(1,INODE)
       LGASH=NGASH+1
       NGASH=LGASH
       BSBAR(1,LGASH)=C2*CARTD(1,INODE)+TEMP1
C
       IF(NTYPE.NE.5)THEN
        TEMP2=C1*CMEAN(2,INODE)
        MGASH=LGASH+1
        NGASH=MGASH
        BSBAR(1,MGASH)=TEMP2-C1*CARTD(2,INODE)
        BSBAR(2,LGASH)=TEMP1-C1*CARTD(1,INODE)
        BSBAR(2,MGASH)=C2*CARTD(2,INODE)+TEMP2
        BSBAR(3,LGASH)=CARTD(2,INODE)
        BSBAR(3,MGASH)=CARTD(1,INODE)
        BSBAR(4,LGASH)=TEMP1-C1*CARTD(1,INODE)
        BSBAR(4,MGASH)=TEMP2-C1*CARTD(2,INODE)
       ENDIF
C
       IF(NTYPE.EQ.3)THEN
        BSBAR(1,LGASH)=BSBAR(1,LGASH)+C1*CMEAN(3,INODE)-
     .                 C1*SHAPE(INODE)/GPCOD(1)
        BSBAR(2,LGASH)=BSBAR(2,LGASH)+C1*CMEAN(3,INODE)-
     .                 C1*SHAPE(INODE)/GPCOD(1)
        BSBAR(4,LGASH)=BSBAR(4,LGASH)+C2*SHAPE(INODE)/GPCOD(1)+
     .                 C1*CMEAN(3,INODE)
       ENDIF           
C
       IF(NTYPE.EQ.4)THEN
        TEMP3=C1*CMEAN(3,INODE)
        NGASH=MGASH+1
        BSBAR(1,NGASH)=TEMP3-C1*CARTD(3,INODE)
        BSBAR(2,NGASH)=TEMP3-C1*CARTD(3,INODE)
        BSBAR(3,NGASH)=0.0D0
        BSBAR(4,NGASH)=C2*CARTD(3,INODE)+TEMP3
        BSBAR(5,LGASH)=CARTD(3,INODE)
        BSBAR(5,MGASH)=0.0D0
        BSBAR(5,NGASH)=CARTD(1,INODE)
        BSBAR(6,LGASH)=0.0D0
        BSBAR(6,MGASH)=CARTD(3,INODE)
        BSBAR(6,NGASH)=CARTD(2,INODE)
       ENDIF
C
      ENDDO
      RETURN
C
      END