      SUBROUTINE CALCPD(UNOMA,STRAT,DEETH,STRRB,STRRC,STRRE,
     .                  RCGTT)
#ifndef restricted
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE PHASE-CHANGE DEFORMATION
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
      INCLUDE 'auxl_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'prob_om.f'
C
      DIMENSION UNOMA(*), STRAT(*), DEETH(*), UNOMX(6)
      DIMENSION RCGTT(*)
C
      IF(NPLATM.EQ.0) RETURN
C
      IF(LARGE.EQ.3) RETURN
C
      ISTAN=2
      IF(LARGE.EQ.0) THEN
       ISTAN=1
      ELSE
       IF(IFREN.EQ.1.OR.IFREN.EQ.4.OR.IFREN.EQ.9.OR.IFREN.EQ.10) ISTAN=1
      ENDIF
C
C**** COMPUTES PHASE-CHANGE DEFORMATION
C
      FTRCX=1.0D0
      FTRBX=1.0D0
C
      IF(ISTAN.EQ.1) THEN                       ! Green-Lagrange strains
       IF(LARGE.NE.0) THEN
        IF(IFREN.EQ.4) THEN
         STRRC=3.0D0*STRRC                      ! volumetric deformation
         STRRB=3.0D0*STRRB
         FTRCX=1.0D0+STRRC
         FTRBX=1.0D0+STRRB
         STRRC=0.5D0*(FTRCX**(2.0D0/3.0D0)-1.0D0)
         STRRB=0.5D0*(FTRBX**(2.0D0/3.0D0)-1.0D0)
        ENDIF
       ENDIF
       DO ISTR1=1,NSTR1
        UNOMX(ISTR1)=UNOMA(ISTR1)
       ENDDO
      ELSE                     ! istan=2         (Almansi strains)
       IF(IFREN.EQ.5.OR.IFREN.EQ.6) THEN
        STRRC=3.0D0*STRRC                       ! volumetric deformation
        STRRB=3.0D0*STRRB
        FTRCX=1.0D0-STRRC
        FTRBX=1.0D0-STRRB
        STRRC=0.5D0*(1.0D0-FTRCX**(2.0D0/3.0D0))
        STRRB=0.5D0*(1.0D0-FTRBX**(2.0D0/3.0D0))
       ENDIF
       DO ISTR1=1,NSTR1           ! transforms to Green-Lagrange strains
        UNOMX(ISTR1)=RCGTT(ISTR1)
       ENDDO
      ENDIF                       ! istan.eq.1
C
      DO ISTR1=1,NSTR1
       STRAT(ISTR1)=STRRB*UNOMX(ISTR1)
      ENDDO
C
C**** COMPUTES INCREMENTAL PHASE-CHANGE DEFORMATION
C
C     Notes:
C
C     For ISTAN=2, DEETH only considers the phase-change function
C     variation of the phase-change deformation, i.e.,
C     dot(E^th) = d STRRB/d fpc * C * dot fpc + STRRB * dot C =
C               = d STRRB/d fpc * C * d fpc/d T * dot T + STRRB * dot C,
C     where DEETH = d STRRB/d fpc * C * d fpc/d T * dot T =
C                 = (STRRB-STRRC) / DTIME * C
C     ( d STRRB/d fpc * dot fpc = (STRRB-STRRC) / DTIME = 
C                               = (STRRB-STRRC) / DTIME * DFPCT, with 
C       DFPCT = d fpc/d T) and the factor DTIME is included in plheat.f
C     This fact is considered in calcbt.f
C
C     Current version: d fpc/d T is assumed to be zero in DEETH. 
C     Besides, the phase-change contribution in DEETH for NITERC=4 is
C     zero as d fpc/d T is assumed to be zero in DEETH
C     (STRRE = d STRRB/d fpc).
C
      IF(ITERME.GT.0) THEN                          ! bidirect. problems
       IF(ITERMP.GT.0) THEN
        IF(NITERC.EQ.1.OR.NITERC.EQ.2.OR.NITERC.EQ.4) THEN
         DO ISTR1=1,NSTR1
          DEETH(ISTR1)=0.0D0                  ! to be improved
c         DEETH(ISTR1)=(STRRB-STRRC)*DFPCT*UNOMX(ISTR1)
          IF (NITERC.EQ.4) THEN
           DEETH(ISTR1)=DEETH(ISTR1)          ! to be improved
c          DEETH(ISTR1)=DEETH(ISTR1)+2.0D0*STRRE*DFPCT*
c    .                    (STRAN(ISTR1)-STRAP(ISTR1)-STRAT(ISTR1))*DTEMG
          ENDIF                ! niterc.eq.4
         ENDDO
        ENDIF                  ! niterc.eq.1.or.niterc.eq.2.or....
       ENDIF                   ! itermp.gt.0
      ENDIF                    ! iterme.gt.0
C
      RETURN
#endif
      END
