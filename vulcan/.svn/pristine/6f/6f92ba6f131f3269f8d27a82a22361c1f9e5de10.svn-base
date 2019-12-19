      SUBROUTINE SUBINC(TEMPN,DTEMP,FPCHA,PWORK,PWORKT,INDEI)
C***********************************************************************
C
C**** THIS ROUTINE SUBINCREMENTS THE TEMPERATURE EFFECT TO BE CONSIDERED
C     IN THE MECHANICAL PROBLEM (INDEI=0; ONLY FOR ITERME=0,2,3,4,5 & 6)
C     AND OBTAINS THE COUPLING TERM DERIVED FROM TEMPERATURE
C     SUBINCREMENTS TO BE CONSIDERED IN THE THERMAL PROBLEM (INDEI=1;
C     ONLY FOR ITERME=2,3,4,5 & 6)
C
C     Note: subinct.f (not implemented yet; Jan/97) should consider the
C           gap effect on the thermal problem for ITERME=1,6
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
      INCLUDE 'inte_om.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION TEMPN(NPOIN,2),     DTEMP(NPOIN),
     .          FPCHA(NFPCH,NPOIN)
C
      DIMENSION PWORK(NPOIN,2),     PWORKT(NPOINT,3)
C
C**** ASSIGNS TEMPN & FPCHA TO LAST CONVERGED TIME STEP VALUES
C
C     TEMPN(1,NPOINC:1,1): temperature
C     FPCHA(1,NNUPT:1,NPOIC): phase-change function
C
C     DTEMP(1,NPOINC): temperature increment
C     FPCHA(1+NNUPT,2*NNUPT:1,NPOIC): phase-change function rate
C
C     TEMPN(1,NPOINC:2,2): initial temperature
C
      IF(INDEI.EQ.0) THEN
       IF(ISSTEP.EQ.1) THEN
        DO IPOINC=1,NPOINC
         TEMPN(IPOINC,1)=TEMPN(IPOINC,1)-DTEMP(IPOINC)
         DTEMP(IPOINC)=DTEMP(IPOINC)/NSSTEP     ! necessary for NITERC>0
         DO IFPCH=1,NNUPT
          FPCHA(IFPCH,IPOINC)=FPCHA(IFPCH,IPOINC)-
     .                        FPCHA(IFPCH+NNUPT,IPOINC)*DTIMET
         ENDDO
        ENDDO
       ENDIF
C
C**** SUBINCREMENTS TEMPN & FPCHA
C
       DO IPOINC=1,NPOINC
        TEMPN(IPOINC,1)=TEMPN(IPOINC,1)+DTEMP(IPOINC)
        DO IFPCH=1,NNUPT
         FPCHA(IFPCH,IPOINC)=FPCHA(IFPCH,IPOINC)+
     .                       FPCHA(IFPCH+NNUPT,IPOINC)*DTIMET/NSSTEP
        ENDDO
       ENDDO
C
      ELSE               ! indei=1
C
C**** OBTAINS THE COUPLING TERM DERIVED FROM TEMPERATURE SUBINCREMENTS
C
       IF(ITERMP.GT.0.OR.ITERMF.GT.0) THEN
        IF(ISSTEP.EQ.1) THEN
         DO IPOINC=1,NPOINC
          PWORKT(IPOINC,1)=PWORK(IPOINC,1)
          IF(NITERC.EQ.2) PWORKT(IPOINC,3)=PWORK(IPOINC,2)
         END DO
        ELSE
         DO IPOINC=1,NPOINC
          PWORKT(IPOINC,1)=PWORKT(IPOINC,1)+PWORK(IPOINC,1)
          IF(NITERC.EQ.2) PWORKT(IPOINC,3)=PWORKT(IPOINC,3)+
     .                                                   PWORK(IPOINC,2)
         END DO
        ENDIF
        IF(ISSTEP.EQ.NSSTEP) THEN        ! to avoid changes in trasmec.f
         DO IPOINC=1,NPOINC
          PWORK(IPOINC,1)=PWORKT(IPOINC,1)/NSSTEP
          IF(NITERC.EQ.2) PWORK(IPOINC,2)=PWORKT(IPOINC,3)/NSSTEP ! ojo
         END DO
        ENDIF
       ENDIF
      ENDIF              ! indei.eq.0
C
      RETURN
      END
