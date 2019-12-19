      SUBROUTINE SETDATD
C***********************************************************************
C
C**** THIS ROUTINE SETS COUPLING DATA (thermal & microstructural prob.)
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'nued_om.f'
C
C**** COUIND
C
      IMICR=1
      IMICO=1      ! coupling index (=0: weak; 1=full)
C
      NPROPM=500   ! microstructural properties; see setdats.f
C
      NNUM4TS=300                             ! change if necessary
c     NNUM4TMS=11*NNUM4TS+46                  ! estimated; see pointes.f
      NNUM4TMS=3*NNUM4TS+6                    ! estimated; see pointes.f
      NHISTM=NNUM4TMS
C
      NNUPM=0           ! number microscopic phase-changes
      NNUPO=0           ! number of other microst. variables
C
      NPLAC=23          ! change if necessary; see nued_om.f
C
      ISTAGGS=0         ! it should be read in couind.f (to be implem.!)
C
      NITERCS=0
C
C***********************************************************************
C
C**** LOGICAL UNITS (defaults)        ! not implemented
C
C***********************************************************************
C
c     LUDATC= 36
c     LURESC=236
C
      RETURN
      END
