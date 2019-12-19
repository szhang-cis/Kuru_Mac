      SUBROUTINE ZEROTE(ELPRE,ELVAR,DISTO,HEADS,TLOAD,RLOAD,REFOR,TEMPN,
     .                  DTEMP,PWORK,PREAS,TGAPS,VNORM,FPCHA,LACTI,PRESF)
C***********************************************************************
C
C**** THIS ROUTINE INITIALISES VARIOUS ARRAYS TO ZERO
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION DISTO(*),           HEADS(*),
     .          TLOAD(*),           RLOAD(*),
     .          REFOR(*)
      DIMENSION ELPRE(NPREV),       ELVAR(NSTAT)
      DIMENSION TEMPN(NPOIN,2),     DTEMP(NPOIN),
     .          PWORK(NPOIN,2)
      DIMENSION PREAS(NPREA,NPOIN), TGAPS(NPOIN), VNORM(NTOTV),
     .          FPCHA(NFPCH,NPOIN), LACTI(NELEM)
      DIMENSION PRESF(NNODE,NDIME,NELEM)
C
      DO 20 INDEX=1,NTOTV*2
      REFOR(INDEX)=0.0D0
   20 TLOAD(INDEX)=0.0D0
C
      DO 30 INDEX=1,NTOTV*NDISO
   30 DISTO(INDEX)=0.0D0
C
      IF(KPORE.NE.0) THEN
       DO 40 INDEX=1,NPOIN*4
   40  HEADS(INDEX)=0.0D0
      ENDIF
C
      IF(ITERME.GE.0) THEN
       DO IPOIN=1,NPOIN
        DO IIINI=1,2
         TEMPN(IPOIN,IIINI)=0.0D0
        ENDDO
        DTEMP(IPOIN)=0.0D0
        DO IFPCH=1,NFPCH
         FPCHA(IFPCH,IPOIN)=0.0D0
        ENDDO
       ENDDO
      ENDIF
C
      IF(ITERME.GT.0) THEN
       DO IPOIN=1,NPOIN
        PWORK(IPOIN,1)=0.0D0
        IF(NITERC.EQ.2.OR.NITERC.EQ.3) PWORK(IPOIN,2)=0.0D0
       ENDDO
      ENDIF
C
      IF(ITERME.GE.0) THEN   ! to print in uni or bidirec. coupled prob.
       DO IPOIN=1,NPOIN
        DO IPREA=1,NPREA
         PREAS(IPREA,IPOIN)=0.0D0
        ENDDO
        TGAPS(IPOIN)=0.0D0
       ENDDO
      ENDIF
C
      DO ITOTV=1,NTOTV
       VNORM(ITOTV)=0.0D0
       RLOAD(ITOTV)=0.0D0
      END DO
C
      IF(NACTI.EQ.1) THEN
       DO IELEM=1,NELEM
        LACTI(IELEM)=0
       ENDDO
      ENDIF
C
      IF(NLDSF.EQ.1) THEN
       DO IELEM=1,NELEM
        DO INODE=1,NNODE
         DO IDIME=1,NDIME
          PRESF(INODE,IDIME,IELEM)=0.0D0
         END DO
        END DO
       END DO
      ENDIF
C
C**** INITIALIZE ELPRE & ELVAR AND WRITE TO DATA BASE 
C
      IF(NMEMOM.EQ.1.OR.INITV.EQ.1) THEN
       DO KPREV=1,NPREV
        ELPRE(KPREV)=0.0D0
       ENDDO
      ENDIF
C
      DO 50 KSTAT=1,NSTAT
   50 ELVAR(KSTAT)=0.0D0
C
c     IF(INITI.EQ.0) RETURN
C
      DO 100 IELEM=1,NELEM
C
       IF(NMEMOM.EQ.1)
     .  CALL DATBAS(ELPRE,    2,    1)
       CALL DATBAS(ELVAR,    3,    1)
C
  100  CONTINUE
C        
      RETURN
      END
