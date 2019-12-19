      SUBROUTINE PLSTRE(VECTG,STRAP,SPTOT,DEEPP,TLAMD,DLAMD,
     .                  DMATX,SIGMA,SGTOT,DESIG,PRESG,
     .                     DS,   DP,
     .                  DTIMX,DAMAX,DAMAY,NAUXI)
c    .                  STRAN,
c    .                  XJACM,XJACI,XJA3M,XJA3I,DETJM,
c    .                  BULKM,DISTM,STRRA,STRR4,
c    .                  RCGTT,RCGTI,TRACC,TRACI,ENEYE)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE STRESS TENSOR
C
C
C     Note:
C
C     IINCO=0 >> Iterative integration
C     IINCO=1 >> Incremental integration
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'auxl_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'prob_om.f'
C
      DIMENSION VECTG(*),       STRAP(*), SPTOT(*), DEEPP(*)
      DIMENSION DMATX(NAUXI,*), SIGMA(*), SGTOT(*), DESIG(*), PRESG(*)
      DIMENSION DS(*),          DP(*)
c     DIMENSION STRAN(*)
c     DIMENSION UNOMA(6)
c     DIMENSION XJACM(NDIME,*), XJACI(NDIME,*),
c    .          RCGTT(*),       RCGTI(*),
c    .          ENEYE(*)
C
      IF(LARGE.EQ.3) RETURN
C
      ISTAN=2
      IF(LARGE.EQ.0) THEN
       ISTAN=1
      ELSE
       IF(IFREN.EQ.1.OR.IFREN.EQ.2.OR.IFREN.EQ.3.OR.
     .    IFREN.EQ.4.OR.IFREN.EQ.5.OR.IFREN.EQ.6.OR.
     .    IFREN.EQ.9.OR.IFREN.EQ.10) ISTAN=1
      ENDIF
C
C**** COMPUTES PLASTIC DEFORMATION AND STRESSES
C
      IF(ISTAN.EQ.1) THEN           ! standard models
C 
C**** COMPUTATION OF THE INCREMENTAL & TOTAL PLASTIC DEFORMATION
C
       IINCO=1                      ! better as input
       DO ISTR1=1,NSTR1
        IF(IINCO.EQ.0) THEN
         DP(ISTR1)=DLAMD*VECTG(ISTR1)*DTIMX
         DEEPP(ISTR1)=DEEPP(ISTR1)+DP(ISTR1)
        ELSE
         DEEPP(ISTR1)=TLAMD*VECTG(ISTR1)*DTIMX
        ENDIF
        SPTOT(ISTR1)=STRAP(ISTR1)+DEEPP(ISTR1)
       ENDDO
C
C**** COMPUTES "PLASTIC STRESS INCREMENT"
C
       DO ISTRS=1,NSTRS
        DS(ISTRS)=0.0D0
       ENDDO
       DO ISTRS=1,NSTRS
        DO JSTRS=1,NSTRS
         DS(ISTRS)=DS(ISTRS)+DMATX(ISTRS,JSTRS)*DEEPP(JSTRS)
        ENDDO
       ENDDO
C
C**** COMPUTES TOTAL & INCREMENTAL STRESSES
C
       DO ISTRS=1,NSTRS
        SGTOT(ISTRS)=(1.0D0-DAMAX)*(SIGMA(ISTRS)-ALFAP*DS(ISTRS))
       ENDDO
      ELSE                          ! istan=2 (non-standard models)
       IF(IFREN.EQ.7) THEN          ! Simo's model
c       CALL RUNEND('ERROR: IFREN=7 NOT IMPLEMENTED')

c       to be revised together with calcst.f !!!!!!!!!!!!!!!!!!!!!!!!

C
C**** COMPUTATION OF THE INCREMENTAL & TOTAL PLASTIC DEFORMATION
C
        IINCO=1                      ! better as input
        DO ISTR1=1,NSTR1
         IF(IINCO.EQ.0) THEN
          DP(ISTR1)=DLAMD*VECTG(ISTR1)*DTIMX
          DEEPP(ISTR1)=DEEPP(ISTR1)+DP(ISTR1)
         ELSE
          DEEPP(ISTR1)=TLAMD*VECTG(ISTR1)*DTIMX
         ENDIF
         SPTOT(ISTR1)=STRAP(ISTR1)+DEEPP(ISTR1)
        ENDDO
C
C**** COMPUTES "PLASTIC STRESS INCREMENT"
C
        DO ISTRS=1,NSTRS
         DS(ISTRS)=0.0D0
        ENDDO
        DO ISTRS=1,NSTRS
         DO JSTRS=1,NSTRS
          DS(ISTRS)=DS(ISTRS)+DMATX(ISTRS,JSTRS)*DEEPP(JSTRS)
         ENDDO
        ENDDO
C
C**** COMPUTES TOTAL & INCREMENTAL STRESSES
C
        DO ISTRS=1,NSTRS
         SGTOT(ISTRS)=(1.0D0-DAMAX)*(SIGMA(ISTRS)-ALFAP*DS(ISTRS))
        ENDDO


        IF(NTYPE.EQ.1) STRR4=0.0D0
       ENDIF            ! ifren.eq.7
       IF(IFREN.EQ.8) THEN          ! idem IFREN=2
        CALL RUNEND('ERROR: IFREN=8 NOT IMPLEMENTED')
        IF(NTYPE.EQ.1) STRR4=0.0D0
       ENDIF            ! ifren.eq.8
      ENDIF             ! istan.eq.1
C
      RETURN
      END
