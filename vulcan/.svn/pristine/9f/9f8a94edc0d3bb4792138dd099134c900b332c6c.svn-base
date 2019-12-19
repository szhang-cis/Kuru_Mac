      SUBROUTINE SSOLID(DVOLU,STRSG)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATE SMOOTHED STRESSES FOR THE 
C     SOLIDIFICATION PROBLEM
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION STRSG(NSTR1,*),DVOLU(*),STSMO(4)
C
      VOLIN=DVOLU(1)+DVOLU(2)+DVOLU(3)+DVOLU(4)
      DO ISTR1=1,NSTR1
       STSMO(ISTR1)=0.0D0
       DO IGAUS=1,NGAUL
        STSMO(ISTR1)=STSMO(ISTR1)+STRSG(ISTR1,IGAUS)*DVOLU(IGAUS)
       ENDDO
       STSMO(ISTR1)=STSMO(ISTR1)/VOLIN
      ENDDO
C
      DO ISTR1=1,NSTR1
       DO IGAUS=1,NGAUL
        STRSG(ISTR1,IGAUS)=STSMO(ISTR1)
       ENDDO
      ENDDO
C
      RETURN
      END
