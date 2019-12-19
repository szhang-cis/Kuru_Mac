      SUBROUTINE TRASTER(DISTOT,TEMPIT,FPCHAT,PWORKT,
     .                    TEMPN, DTEMP, FPCHA, PWORK)
C***********************************************************************
C
C**** THIS ROUTINE TRANSFER NODAL THERMAL VARIABLES TO MECHANICAL
C     PROBLEM
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
      INCLUDE 'prob_om.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
C
      DIMENSION DISTOT(NTOTVT,3),     TEMPIT(NPOINT,2),
     .          FPCHAT(NFPCH,NPOINT), PWORKT(NPOINT,3),
     .          TEMPN(NPOIN,2),       DTEMP(NPOIN),
     .          FPCHA(NFPCH,NPOIN),   PWORK(NPOIN,2)
C
      DO IPOINC=1,NPOINC
       TEMPAUX=TEMPN(IPOINC,1)
       TEMPN(IPOINC,1)=DISTOT(IPOINC,1)
       TEMPN(IPOINC,2)=TEMPIT(IPOINC,1)            ! initial temperature
       IF(NSSTEPT.EQ.1) THEN
        DTEMP(IPOINC)=DISTOT(IPOINC,2)*DTIMET
       ELSE
        DTEMP(IPOINC)=TEMPN(IPOINC,1)-TEMPAUX
       ENDIF
       DO IFPCH=1,NFPCH
        FPCHA(IFPCH,IPOINC)=FPCHAT(IFPCH,IPOINC)
       ENDDO
      ENDDO
C
      IF(ITERMP.GT.0) THEN
       IF(NITERC.EQ.3) THEN    ! adiabatic staggered scheme
        DO IPOINC=1,NPOINC
         PWORK(IPOINC,2)=PWORKT(IPOINC,3)
        END DO
       ENDIF
      ENDIF
C
      RETURN
      END
