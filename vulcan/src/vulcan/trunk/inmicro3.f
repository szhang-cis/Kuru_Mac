      SUBROUTINE INMICRO3(ALPHAM,INUPC,IPLAT)
C***********************************************************************
C
C**** THIS ROUTINE INITIALITES MICROSTRUCTURAL PARAMETERS FOR MODEL 3
C
C***********************************************************************
C
C     Index of variables
C
C     ALPHAM= array of microstructural (microscopical) variables
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** COUPLING VARIABLES (thermal-microstructural)
C
      INCLUDE 'nued_om.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION ALPHAM(NHISTM)
C
C**** TRANSFERS "VPLAT" TO MICROSCOPICAL VARIABLES
C
      SICON=VPLAT(IPLAT,6)
      SBCON=VPLAT(IPLAT,7)
      SICONE=VPLAT(IPLAT,8)
      NLICO=INT(VPLAT(IPLAT,9))
      PARCO=VPLAT(IPLAT,10)
C
      IMNUC=INT(VPLAT(IPLAT,11))
      IMGRO=INT(VPLAT(IPLAT,12))
      IMNUCE=INT(VPLAT(IPLAT,13))
      IMGROE=INT(VPLAT(IPLAT,14))
C
      IF(IMNUC.EQ.1.OR.IMNUC.EQ.2) THEN
       TEMAV=VPLAT(IPLAT,14+1)
       TEMDE=VPLAT(IPLAT,14+2)
       GRDEM=VPLAT(IPLAT,14+3)
       IVN=3
      ENDIF
C
      IF(IMGRO.EQ.1.OR.IMGRO.EQ.2.OR.IMGRO.EQ.3.OR.IMGRO.EQ.4.OR.
     .   IMGRO.EQ.5) THEN
       NRAUN=INT(VPLAT(IPLAT,14+IVN+1))
       DIFLI=VPLAT(IPLAT,14+IVN+2)
       HENER=VPLAT(IPLAT,14+IVN+3)
       IVG=3
      ENDIF
C
      IF(IMNUCE.EQ.1.OR.IMNUCE.EQ.2) THEN
       TEMAVE=VPLAT(IPLAT,14+IVN+IVG+1)
       TEMDEE=VPLAT(IPLAT,14+IVN+IVG+2)
       GRDEME=VPLAT(IPLAT,14+IVN+IVG+3)
       IVNE=3
      ENDIF
C
      IF(IMGROE.EQ.1.OR.IMGROE.EQ.2.OR.IMGROE.EQ.3.OR.
     .   IMGROE.EQ.4.OR.IMGROE.EQ.5.OR.IMGROE.EQ.6) THEN
       NRAUNE=INT(VPLAT(IPLAT,14+IVN+IVG+IVNE+1))
       HENERE=VPLAT(IPLAT,14+IVN+IVG+IVNE+2)
       IVGE=2
      ENDIF
C
      IFPCDT= INT(VPLAT(IPLAT,14+IVN+IVG+IVNE+IVGE+1))
      IAFLOJ= INT(VPLAT(IPLAT,14+IVN+IVG+IVNE+IVGE+2))
      NSUGUEM=INT(VPLAT(IPLAT,14+IVN+IVG+IVNE+IVGE+3))
      EXTGUE=     VPLAT(IPLAT,14+IVN+IVG+IVNE+IVGE+4)
      IVFPC=4
C
      IV=                        ! number of VPLAT defined in input data
     .   14+IVN+IVG+IVNE+IVGE+IVFPC
C
C**** INITIALITES MICROSTRUCTURAL PARAMETERS ("ALPHAM" & OTHERS)
C
      ALPHAM(1)=1.0D0     ! liquid state is assumed as initial condition
      DO IN=2,14
       ALPHAM(IN)=0.0D0
      ENDDO
C
C**** INCREMENTS "ALPHAM" INDEX
C
      INUPC=INUPC+14                       ! less than 19; see pointes.f
C
      RETURN
      END
