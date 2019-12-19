      SUBROUTINE ASSPRO(PROPS,PROPSS)
C***********************************************************************
C
C**** THIS ROUTINE ASSIGNS THE MICROSTRUCTURAL PROPERTIES TO MECHANICAL
C     PROPERTIES
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nued_om.f'                      ! thermal-microstructural
      INCLUDE 'nuee_om.f'                   ! mechanical-microstructural
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'auxl_om.f'
      INCLUDE 'prob_om.f'
C
C**** MICROSCOPICAL VARIABLES
C
      INCLUDE 'auxl_oms.f'
      INCLUDE 'prob_oms.f'
C
      DIMENSION PROPS(NPROP,*), PROPSS(NPROPS,*)
C
C**** ASSIGNS MATERIALS PROPERTIES
C
      NPROPX=NPROP-NPROPM           ! number of microstructural prop.
      DO IMATS=1,NMATS
       DO IPROP=1,NPROPX
        PROPS(IPROP+NPROPM,IMATS)=PROPSS(IPROP,IMATS)
       ENDDO
      ENDDO
C
      RETURN
      END
