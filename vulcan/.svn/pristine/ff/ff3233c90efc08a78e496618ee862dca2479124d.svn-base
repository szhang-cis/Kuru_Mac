      SUBROUTINE CALCLD(UNOMA,STRAT,STRAN,STRAP,RCGTT,
     .                  XJACM,XJACI,XJA3M,XJA3I,
     .                  SHRIN)
#ifndef restricted
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE LIQUID & SHRINKAGE DEFORMATIONS
C
C     In the liquid, the plastic deformation is assumed to be equals
C     to the deviatoric total deformation.
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
      DIMENSION UNOMA(*),       STRAT(*)
      DIMENSION STRAN(*),       STRAP(*),       RCGTT(*)
      DIMENSION XJACM(NDIME,*), XJACI(NDIME,*)
      DIMENSION STRANL(6),      STRAPL(6)


c     return                       ! to be revised (July 1999)


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
C**** COMPUTES THE LIQUID DEFORMATION
C
      IF(ISTAN.EQ.1) THEN                       ! Green-Lagrange strains
       DEFVO=STRAN(1)+STRAN(2)+STRAN(4)         ! volumetric strain
       DO ISTRS=1,NSTRS
        STRAP(ISTRS)=STRAN(ISTRS)-1.0D0/3.0D0*DEFVO
        IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6)
     .                                         STRAP(ISTRS)=STRAN(ISTRS)
       ENDDO
      ELSE                     ! istan=2         (Almansi strains)
       IF(IFREN.EQ.2.OR.IFREN.EQ.3.OR.
     .    IFREN.EQ.5.OR.IFREN.EQ.6.OR.IFREN.EQ.8) THEN
        DETJX=1.0
        CALL PUFOBA(STRAN,STRANL,XJACI,XJA3I,DETJX,    2)
        DEFVO=STRANL(1)+STRANL(2)+STRANL(4)     ! volumetric strain
        DO ISTRS=1,NSTRS
         STRAPL(ISTRS)=STRANL(ISTRS)-1.0D0/3.0D0*DEFVO
         IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6)
     .                                       STRAPL(ISTRS)=STRANL(ISTRS)
        ENDDO
        CALL PUFOBA(STRAPL,STRAP,XJACM,XJA3M,DETJX,    2)
       ENDIF
       IF(IFREN.EQ.7) THEN
        CALL RUNEND('ERROR: MODEL 7 NOT IMPLEMENTED IN CALCLD')
       ENDIF
      ENDIF                    ! istan.eq.1
C
C**** COMPUTES THE SHRINKAGE DEFORMATION
C
      IF(ISHRI.EQ.1) THEN
       IF(ISTAN.EQ.1) THEN                      ! Green-Lagrange strains
        SHRIA=0.0D0
        DO ISTR1=1,NSTR1
         SHRIA=SHRIA+(STRAN(ISTR1)-STRAT(ISTR1))*UNOMA(ISTR1)
        ENDDO
        SHRIN=0.0D0                         ! reversibility
        IF(SHRIA.GT.0.0) SHRIN=SHRIA
       ELSE                    ! istan=2         (Almansi strains)
        IF(IFREN.EQ.2.OR.IFREN.EQ.3.OR.
     .     IFREN.EQ.5.OR.IFREN.EQ.6.OR.IFREN.EQ.8) THEN
         CALL RUNEND('ERROR: MODEL 2,3,5,6,8 NOT IMPLEMENTED IN CALCLD')
        ENDIF
        IF(IFREN.EQ.7) THEN
         CALL RUNEND('ERROR: MODEL 7 NOT IMPLEMENTED IN CALCLD')
        ENDIF
       ENDIF                   ! istan.eq.1
      ENDIF                    ! ishri.eq.1
C
      RETURN
#endif
      END
