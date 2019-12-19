      SUBROUTINE INCS01(EHIST,PROPS,STRAN,TEMPC)
C***********************************************************************
C
C*** THIS ROUTINE INCREMENTS NON-TENSIONAL STRAINS AT GAUSS POINTS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      COMMON/INCSTNA/IFUNC
      COMMON/INCSTNB/HFUNC(5)
C
      DIMENSION PROPS(*), EHIST(NHIST,*), STRAN(NSTR1,NGAUS), TEMPC(*)
C
C***DEAL WITH TEMPERATURE
C
      IF(KTEMP.EQ.1) THEN
        ALPHA=PROPS(43)
        FACTT=TEMPC(1)+TEMPC(2)*DSIN(TEMPC(3)+TEMPC(4)*TTIME)
        IF(ITIME*ISTEP.EQ.1) THEN
          FACTP=0.0
        ELSE
          FACTP=TEMPC(1)+TEMPC(2)*DSIN(TEMPC(3)+TEMPC(4)*(TTIME-DTIME))
        ENDIF
        DELTA=(FACTT-FACTP)*ALPHA
        IF(DELTA.EQ.0.0) RETURN
C
C       increment strains
C
        DO IGAUS=1,NGAUL
                         STRAN(1,IGAUS)=STRAN(1,IGAUS)+DELTA
                         STRAN(2,IGAUS)=STRAN(2,IGAUS)+DELTA
          IF(NSTR1.GE.4) STRAN(4,IGAUS)=STRAN(4,IGAUS)+DELTA
        ENDDO
C
      ENDIF
C
C***DEAL WITH PRESCRIBED STRAINS
C
      ITYPE=INT(HFUNC(1))
      IF(IFUNC.EQ.0.OR.ITYPE.EQ.0) RETURN
C
      IF(IFUNC.GT.0) THEN            ! SAME CURVE FOR ALL POINTS    
C
C       current volumetric strains
C
        FACTT=0.0
             IF(ITYPE.EQ.1) THEN
                            CALL FUNLO1(FACTT,HFUNC,TTIME)
        ELSE IF(ITYPE.EQ.2) THEN
                            CALL FUNLO2(FACTT,HFUNC,TTIME)
        ELSE IF(ITYPE.EQ.3) THEN
                            CALL FUNLO3(FACTT,HFUNC,TTIME)
        ELSE IF(ITYPE.EQ.4) THEN 
                            CALL FUNLO4(FACTT,HFUNC,TTIME)
        ENDIF
C
C       previous volumetric strains
C
        FACTP=0.0
             IF(ITYPE.EQ.1) THEN
                            CALL FUNLO1(FACTP,HFUNC,TTIME-DTIME)
        ELSE IF(ITYPE.EQ.2) THEN
                            CALL FUNLO2(FACTP,HFUNC,TTIME-DTIME)
        ELSE IF(ITYPE.EQ.3) THEN
                            CALL FUNLO3(FACTP,HFUNC,TTIME-DTIME)
        ELSE IF(ITYPE.EQ.4) THEN 
                            CALL FUNLO4(FACTP,HFUNC,TTIME-DTIME)
        ENDIF
C
        DELTA=FACTT-FACTP
C
        IF(DELTA.EQ.0.0) RETURN
C
C       increment strains
C
        DO 10 IGAUS=1,NGAUL
C
                       STRAN(1,IGAUS)=STRAN(1,IGAUS)+DELTA
                       STRAN(2,IGAUS)=STRAN(2,IGAUS)+DELTA
        IF(NSTR1.GE.4) STRAN(4,IGAUS)=STRAN(4,IGAUS)+DELTA
C
   10   CONTINUE
C
      ELSE                          ! DIFFERENT CURVE FOR EACH POINT
C
        DO 20 IGAUS=1,NGAUL
C
        TINIT=EHIST(19,IGAUS)         ! Initial time for this point
        IF(TINIT.EQ.0.0D+00) GO TO 20
C
        RTIME=TTIME-TINIT             ! 'Corrected' time for this point
C
C       current volumetric strains
C
        FACTT=0.0
             IF(ITYPE.EQ.1) THEN
                            CALL FUNLO1(FACTT,HFUNC,RTIME)
        ELSE IF(ITYPE.EQ.2) THEN
                            CALL FUNLO2(FACTT,HFUNC,RTIME)
        ELSE IF(ITYPE.EQ.3) THEN
                            CALL FUNLO3(FACTT,HFUNC,RTIME)
        ELSE IF(ITYPE.EQ.4) THEN 
                            CALL FUNLO4(FACTT,HFUNC,RTIME)
        ENDIF
C
C       previous volumetric strains
C
        FACTP=0.0
             IF(ITYPE.EQ.1) THEN
                            CALL FUNLO1(FACTP,HFUNC,RTIME-DTIME)
        ELSE IF(ITYPE.EQ.2) THEN
                            CALL FUNLO2(FACTP,HFUNC,RTIME-DTIME)
        ELSE IF(ITYPE.EQ.3) THEN
                            CALL FUNLO3(FACTP,HFUNC,RTIME-DTIME)
        ELSE IF(ITYPE.EQ.4) THEN 
                            CALL FUNLO4(FACTP,HFUNC,RTIME-DTIME)
        ENDIF
C
        DELTA=FACTT-FACTP
C
C       increment strains
C
        if(abs(rtime-dtime).lt.1.0e-08) 
     .                        print 900, igaus,ielem,tinit
  900   format(2x,'El punto',i5,'  del elemento',i6,
     .            ' rompio en ',f15.6)
                       STRAN(1,IGAUS)=STRAN(1,IGAUS)+DELTA
                       STRAN(2,IGAUS)=STRAN(2,IGAUS)+DELTA
        IF(NSTR1.GE.4) STRAN(4,IGAUS)=STRAN(4,IGAUS)+DELTA
C
   20   CONTINUE
C
      ENDIF
C
      RETURN
      END
