      SUBROUTINE HOURB4(BEXI,BETA,CARTD,DAREA,ELCOD,
     .                  NDIME,NEVAB,NNODE,NSTRS)
C***********************************************************************
C
C****THIS ROUTINE COMPURES THE "HOURGLASS" B-MATRICES FOR THE 
C       4-NODE ELEMENT
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION CARTD(NDIME,*), ELCOD(NDIME,*)  
      DIMENSION BEXI(NSTRS,*),  BETA(NSTRS,*)
      DIMENSION GAMMA(4),  EXI(4),    ETA(4),   H(4)
      DATA EXI/-1.0, 1.0, 1.0,-1.0/
      DATA ETA/-1.0,-1.0, 1.0, 1.0/
      DATA   H/ 1.0,-1.0, 1.0,-1.0/
C
C     DEFINE INNER PRODUCTS
C
      HX=0.0
      HY=0.0
      EXIX=0.0
      EXIY=0.0
      ETAX=0.0
      ETAY=0.0
      DO I=1,NNODE
        HX  =  HX+  H(I)*ELCOD(1,I)
        HY  =  HY+  H(I)*ELCOD(2,I)
        EXIX=EXIX+EXI(I)*ELCOD(1,I)
        EXIY=EXIY+EXI(I)*ELCOD(2,I)
        ETAX=ETAX+ETA(I)*ELCOD(1,I)
        ETAY=ETAY+ETA(I)*ELCOD(2,I)
      ENDDO
C
C     DEFINE GAMMA
C
      DO I=1,NNODE
        GAMMA(I)=H(I)-HX*CARTD(1,I)-HY*CARTD(2,I)
      ENDDO
C
C     DEFINE B MATRIX
C
      COE13=1.0/3.0
      COE23=2.0/3.0
      IDEVF=0
      IF(IDEVF.EQ.1) THEN
        COE13=0.0
        COE23=1.0
      ENDIF
      DO I=1,NNODE
        J =(I-1)*2+1
        J1=J+1
        BEXI(1,J )=COE23*(-EXIY)*GAMMA(I)
        BEXI(2,J )=COE13*( EXIY)*GAMMA(I)
        BEXI(3,J )=      ( EXIX)*GAMMA(I)
        BEXI(1,J1)=COE13*(-EXIX)*GAMMA(I)
        BEXI(2,J1)=COE23*( EXIX)*GAMMA(I)
        BEXI(3,J1)=      (-EXIY)*GAMMA(I)
        BETA(1,J )=COE23*( ETAY)*GAMMA(I)
        BETA(2,J )=COE13*(-ETAY)*GAMMA(I)
        BETA(3,J )=      (-ETAX)*GAMMA(I)
        BETA(1,J1)=COE13*( ETAX)*GAMMA(I)
        BETA(2,J1)=COE23*(-ETAX)*GAMMA(I)
        BETA(3,J1)=      ( ETAY)*GAMMA(I)
        IF(NSTRS.EQ.4) THEN
          BEXI(4,J )=COE13*( EXIY)*GAMMA(I)
          BEXI(4,J1)=COE13*(-EXIX)*GAMMA(I)
          BETA(4,J )=COE13*(-ETAY)*GAMMA(I)
          BETA(4,J1)=COE13*( ETAX)*GAMMA(I)
        ENDIF
      ENDDO
      COEFO=0.25/DAREA
      DO ISTRE=1,NSTRS
        DO IEVAB=1,NEVAB
          BEXI(ISTRE,IEVAB)=COEFO*BEXI(ISTRE,IEVAB)
          BETA(ISTRE,IEVAB)=COEFO*BETA(ISTRE,IEVAB)
        ENDDO
      ENDDO
      RETURN
      END
