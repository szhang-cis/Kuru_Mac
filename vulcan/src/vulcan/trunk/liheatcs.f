      SUBROUTINE LIHEATCS(DSTRAS,TSTRAS,DMATXS,CONVES,NDIMES,NSTR1S,
     .                    NNODES,NDOFCS,CARTDS,VELCMS,ELDISS,KDYNAS)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE CONVECTIVE TERM
C     (TEMPERATURE GRADIENT * DENSITY * CAPACITY * VELOCITY)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'auxl_oms.f'
      INCLUDE 'inte_oms.f'
C
      DIMENSION DSTRAS(*),        TSTRAS(*),
     .          DMATXS(NSTR1S,*)
      DIMENSION CARTDS(NDIMES,*), ELDISS(NDOFCS,*),
     .          VELCMS(*)
C
C**** CALCULATE THE TOTAL TEMPERATURE GRADIENT IN TSTRAS
C
      DO IDIMES=1,NDIMES
       TSTRAS(IDIMES)=0.0
       IF(KDYNAS.EQ.1) THEN
        DO IEVABS=1,NNODES
         TSTRAS(IDIMES)=TSTRAS(IDIMES)+CARTDS(IDIMES,IEVABS)*
     .                  ELDISS(NDOFCS,IEVABS)
        END DO
       ENDIF
      END DO
C
C**** CALCULATE THE INCREMENTAL TEMPERATURE GRADIENT IN DSTRAT
C
      DO IDIMES=1,NDIMES
       DSTRAS(IDIMES)=0.0
       DO IEVABS=1,NNODES
        DSTRAS(IDIMES)=DSTRAS(IDIMES)+CARTDS(IDIMES,IEVABS)*
     .                 VELCMS(IEVABS)*DTIMES
       END DO
      END DO
C
C**** CALCULATE THE CONVECTIVE TERM
C
      CONVES=0.0
      DO IDIMES=1,NDIMES
       CONVES=CONVES+DMATXS(IDIMES,1)*TSTRAS(IDIMES)
      ENDDO
C
      IF(KDYNAS.EQ.1) THEN     ! transient
       IF(KINTES.EQ.1) THEN    ! Euler method
        CONVES=0.0
        DO IDIMES=1,NDIMES
         CONVES=CONVES+DMATXS(IDIMES,1)*
     .         (TALFAS*TSTRAS(IDIMES)+
     .         (1.0D+00-TALFAS)*(TSTRAS(IDIMES)-DSTRAS(IDIMES)))
        ENDDO
       ENDIF
      ENDIF
C
      RETURN
      END
