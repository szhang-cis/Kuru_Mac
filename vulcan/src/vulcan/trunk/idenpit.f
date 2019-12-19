      SUBROUTINE IDENPIT(MULPHT,MULPTT,MULPDT,NPHCTT,NPHCDT,NISTTT,
     .                   NISDTT,ELDIST,PROPST,VELCMT,ISINRT,TEMTLT,
     .                   TEMDLT)
C***********************************************************************
C
C**** THIS ROUTINE IDENTIFIES UNI OR MULTIPHASE ELEMENTS
C
C***********************************************************************
C
C     NPLAT=NUMBER OF DIFFERENT PLATEAUX
C     TEINF=SOLIDUS TEMPERATURE
C     TESUP=LIQUIDUS TEMPERATURE
C
C***********************************************************************
C
C     MULPHT=INDEX TO IDENTIFY MULTIPHASE ELEMENTS IN TIME t AND t+dt
C           =1    UNIPHASE ELEMENT
C           =2    MULTIPHASE ELEMENT
C
C     MULPTT=INDEX TO IDENTIFY MULTIPHASE ELEMENTS IN TIME t
C           =1    UNIPHASE ELEMENT
C           =2    BI OR MULTIPHASE ELEMENT
C
C     MULPDT=INDEX TO IDENTIFY MULTIPHASE ELEMENTS IN TIME t+dt
C           =1    UNIPHASE ELEMENT
C           =2    BI OR MULTIPHASE ELEMENT
C
C     NPHCTT=INDEX TO IDENTIFY THE INTEGRATION FORM OF BIPHASE 
C            ELEMENTS IN TIME t
C
C     NPHCDT=INDEX TO IDENTIFY THE INTEGRATION FORM OF BIPHASE 
C            ELEMENTS IN TIME t+dt
C
C     NISOTT=INDEX TO IDENTIFY THE PHASE CHANGE CHARACTERISTIC
C           =1    ISOTHERMAL PHASE CHANGE      >>   BIPHASE ELEMENTS
C           =2    NON-ISOTHERMAL PHASE CHANGE  >>   MULTIPHASE ELEMENTS
C            (NISTTT IN TIME t & NISDTT IN TIME t+dt)
C
C     ISINRT=INDEX USED IN ELM005
C           =1  CONTRIBUTION TO THE CAPACITY MATRIX
C           =2  CONTRIBUTION TO THE "THERMAL FORCES" VECTOR
C
C     TEMPLT=PHASE CHANGE TEMPERATURE IN ISOTHERMAL PHASE CHANGE
C            (TEMTLT IN TIME t & TEMDLT IN TIME t+dt)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nued_om.f'   ! thermal-microstructural
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION ELDIST(NDOFCT,*), PROPST(*),
     .          VELCMT(*)
      DIMENSION TESUP(5),         TEINF(5),
     .          NSOTE(5),         MSOTH(5),
     .          MSOTT(5),         MSOTD(5)
C
      MULPHT=1
      IF(NPLAT.EQ.0) THEN
       MULPTT=1
       MULPDT=1
       RETURN
      ENDIF
C
C**** ESTABLISHES SOME PARAMETERS FOR MICROSTRUCTURAL PROBLEMS
C
      IF(IMICR.EQ.1) THEN    ! microscopical approach
       MULPTT=1              ! the subdomain technique is not available
       MULPDT=1
C
       NPHCTT=1              ! just to assign something
       NISTTT=1
       TEMTLT=0.0
       NPHCDT=1
       NISDTT=1
       TEMDLT=0.0
       RETURN
      ENDIF         ! imicr.eq.1
C
C**** CALCULATE HOW THE INTERNAL ENERGY JUMP IS
C
      DO IPLAT=1,NPLAT
       TEINF(IPLAT)=VPLAT(IPLAT,1)
       TESUP(IPLAT)=VPLAT(IPLAT,2)
      ENDDO
C
      TEMTOL=1.0D-10
      DO IPLAT=1,NPLAT
       NSOTE(IPLAT)=1
       DELTT=TESUP(IPLAT)-TEINF(IPLAT)
       IF(DELTT.LT.0.0) 
     .  CALL RUNENDT('IDENPIT: T_SOLIDUS > T_LIQUIDUS    ')
       DELTT=DABS(DELTT)
       IF(DELTT.GT.TEMTOL) NSOTE(IPLAT)=2
      ENDDO
C
      MULPTT=0
      MULPDT=0
      NISTTT=0
      NISDTT=0
      TEMTLT=0.0
      TEMDLT=0.0
C
      DO IPLAT=1,NPLAT
       MSOTH(IPLAT)=1
       MSOTT(IPLAT)=0
       MSOTD(IPLAT)=0
      ENDDO
C
C**** EVALUATE HOW THE ELEMENT IS IN TIME t
C
C     MULPTT=1   UNIPHASE ELEMENT
C     MULPTT=2   BI OR MULTIPHASE ELEMENT
C
      DO IPLAT=1,NPLAT
       NSOLD=0
       NMUSH=0
       NLIQD=0
       DO INODLT=1,NNODLT
        ELDIF=ELDIST(1,INODLT)-VELCMT(INODLT)*DTIMET
C
        IF(ELDIF.GT.TESUP(IPLAT))       NLIQD=NLIQD+1
        IF((ELDIF.LE.TESUP(IPLAT)).AND.
     .     (ELDIF.GT.TEINF(IPLAT)))     NMUSH=NMUSH+1
        IF(ELDIF.LE.TEINF(IPLAT))       NSOLD=NSOLD+1
       END DO
C
       IF(NSOLD.NE.0) MSOTT(IPLAT)=MSOTT(IPLAT)+1
       IF(NMUSH.NE.0) MSOTT(IPLAT)=MSOTT(IPLAT)+1
       IF(NLIQD.NE.0) MSOTT(IPLAT)=MSOTT(IPLAT)+1
      ENDDO
C
      IBOTO=0
      DO IPLAT=1,NPLAT
       IF(MSOTT(IPLAT).GT.1) THEN
        IBOTO=IBOTO+1
        IBOTX=IPLAT
       ENDIF
      ENDDO
      IF(IBOTO.EQ.0) MULPTT=1
      IF(IBOTO.EQ.1) THEN
       MULPTT=2
       NISTTT=NSOTE(IBOTX)
       TEMTLT=TESUP(IBOTX)
      ENDIF
C
C**** WHEN IBOTO=2 THE ELEMENT CONTAINS TWO OR MORE DIFFERENT PLATEAUXS
C     AT THE SAME TIME t
C     >> WARNING: ADOPT ANOTHER DISCRETIZATION
C
      IF(IBOTO.GE.2) THEN
       CALL RUNMENT('IDENPIT: WARNING WITH IBOTO GE 2')
       MULPTT=2
       NISTTT=2      ! is treated as a non-isothermal phase-change
      ENDIF
C
C**** EVALUATE HOW THE ELEMENT IS IN TIME t + dt
C
C     MULPDT=1   UNIPHASE ELEMENT
C     MULPDT=2   BI OR MULTIPHASE ELEMENT
C
      DO IPLAT=1,NPLAT
       NSOLD=0
       NMUSH=0
       NLIQD=0
       DO INODLT=1,NNODLT
        ELDIF=ELDIST(1,INODLT)
C
        IF(ELDIF.GT.TESUP(IPLAT))       NLIQD=NLIQD+1
        IF((ELDIF.LE.TESUP(IPLAT)).AND.
     .     (ELDIF.GT.TEINF(IPLAT)))     NMUSH=NMUSH+1
        IF(ELDIF.LE.TEINF(IPLAT))       NSOLD=NSOLD+1
       END DO
C
       IF(NSOLD.NE.0) MSOTD(IPLAT)=MSOTD(IPLAT)+1
       IF(NMUSH.NE.0) MSOTD(IPLAT)=MSOTD(IPLAT)+1
       IF(NLIQD.NE.0) MSOTD(IPLAT)=MSOTD(IPLAT)+1
      ENDDO
C
      IBOTO=0
      DO IPLAT=1,NPLAT
       IF(MSOTD(IPLAT).GT.1) THEN
        IBOTO=IBOTO+1
        IBOTX=IPLAT
       ENDIF
      ENDDO
      IF(IBOTO.EQ.0) MULPDT=1
      IF(IBOTO.EQ.1) THEN
       MULPDT=2
       NISDTT=NSOTE(IBOTX)
       TEMDLT=TESUP(IBOTX)
      ENDIF
C
C**** WHEN IBOTO=2 THE ELEMENT CONTAINS TWO OR MORE DIFFERENT PLATEAUXS
C     AT THE SAME TIME t+dt
C     >> WARNING: ADOPT ANOTHER DISCRETIZATION
C
      IF(IBOTO.GE.2) THEN
       CALL RUNMENT('IDENPIT: WARNING WITH IBOTO GE 2')
       MULPDT=2
       NISDTT=2      ! is treated as a non-isothermal phase-change
      ENDIF
C
      IF(LINPC.EQ.1) THEN                  ! deals with linearized ph-ch
       MULPDT=MULPTT
       NPHCDT=NPHCTT
       NISDTT=NISTTT
       TEMDLT=TEMTLT
      ENDIF
C
      IF(ISINRT.EQ.1) THEN
C
C     ISINRT=1  CONTRIBUTION TO THE CAPACITY MATRIX
C
C**** MULPHT=1 UNIPHASE ELEMENT    MULPHT=2 BI OR MULTIPHASE ELEMENT
C
       IF(MULPTT.GT.1.OR.MULPDT.GT.1) MULPHT=2
C
       RETURN
      ENDIF
C
C**** ISINRT=2 !!!!
C     ISINRT=2  CONTRIBUTION TO THE "THERMAL FORCES" VECTOR
C
      IF(MULPTT.EQ.1) GO TO 11
      IF(NISTTT.EQ.2) GO TO 11
C
C**** IDENTIFY INTEGRATION ELEMENT FOR THE INTERNAL ENERGY FUNCTION IN
C     TIME t FOR THE CASE NISTTT=1 (BIPHASE ELEMENT)
C
      NPHCTT=0
C
C**** ADDITIONAL CONTROL
C
      IF(NRULET.GT.5)
     . CALL RUNENDT('ERROR:ISOTH. PHASE-CHANGE WITH IMPROVED INT. RULE')
C
      GOTO(1,2,3), NDIMET
C
    1 CONTINUE
C
C**** 2 NODED ELEMENT
C
      IF(NNODLT.EQ.2) THEN
       NPHCTT=1
       GO TO 11
      ENDIF
C
C**** 3 NODED ELEMENT
C
      IF(NNODLT.EQ.3) THEN
       CALL RUNENDT('IDENPIT:ERROR DETECTED WITH NDIME=1')
c
c      not implemented yet
c
c      ELDI1=ELDIST(1,1)-VELCMT(1)*DTIMET
c      ELDI2=ELDIST(1,2)-VELCMT(2)*DTIMET
c      ELDI3=ELDIST(1,3)-VELCMT(3)*DTIMET
c      IF(ELDI1.GT.ELDI3.AND.ELDI2.GT.ELDI3) NPHCTT=7
c      IF(ELDI1.LT.ELDI3.AND.ELDI2.LT.ELDI3) NPHCTT=7
c      NPHCTT=6
       GO TO 11
      ENDIF
C
    2 CONTINUE
C
C**** 3 NODED ELEMENT
C
      IF(NNODLT.EQ.3) THEN
       ELDI1=ELDIST(1,1)-VELCMT(1)*DTIMET
       ELDI2=ELDIST(1,2)-VELCMT(2)*DTIMET
       ELDI3=ELDIST(1,3)-VELCMT(3)*DTIMET
C
       CALL COIN02T(ELDI2,ELDI1,TEMTLT,RCOOR,ISOL1)
       CALL COIN02T(ELDI3,ELDI1,TEMTLT,RCOOR,ISOL2)
       CALL COIN02T(ELDI3,ELDI2,TEMTLT,RCOOR,ISOL3)
C
       if(isol1.gt.0.and.isol2.eq.0.and.isol3.gt.0) nphctt=7
       if(isol1.gt.0.and.isol2.gt.0.and.isol3.eq.0) nphctt=8
       if(isol1.eq.0.and.isol2.gt.0.and.isol3.gt.0) nphctt=9
C
       if(nphctt.eq.0) mulptt=1      ! to avoid unnecessary subdivisions
C
       GO TO 11
      ENDIF
C
C**** 4 NODED ELEMENT
C
      IF(NNODLT.EQ.4) THEN
       ELDI1=ELDIST(1,1)-VELCMT(1)*DTIMET
       ELDI2=ELDIST(1,2)-VELCMT(2)*DTIMET
       ELDI3=ELDIST(1,3)-VELCMT(3)*DTIMET
       ELDI4=ELDIST(1,4)-VELCMT(4)*DTIMET
C
       CALL COIN01T(ELDI2,ELDI1,TEMTLT,RCOOR,ISOL1)
       CALL COIN01T(ELDI3,ELDI2,TEMTLT,RCOOR,ISOL2)
       CALL COIN01T(ELDI4,ELDI3,TEMTLT,RCOOR,ISOL3)
       CALL COIN01T(ELDI1,ELDI4,TEMTLT,RCOOR,ISOL4)
C
       if(isol1.gt.0.and.isol2.eq.0.and.isol3.eq.0.and.isol4.gt.0) 
     .                                                         nphctt=1
       if(isol1.gt.0.and.isol2.gt.0.and.isol3.eq.0.and.isol4.eq.0) 
     .                                                         nphctt=2
       if(isol1.eq.0.and.isol2.gt.0.and.isol3.gt.0.and.isol4.eq.0) 
     .                                                         nphctt=3
       if(isol1.eq.0.and.isol2.eq.0.and.isol3.gt.0.and.isol4.gt.0) 
     .                                                         nphctt=4
       if(isol1.gt.0.and.isol2.eq.0.and.isol3.gt.0.and.isol4.eq.0) 
     .                                                         nphctt=5
       if(isol1.eq.0.and.isol2.gt.0.and.isol3.eq.0.and.isol4.gt.0) 
     .                                                         nphctt=6
       if(isol1.gt.0.and.isol2.gt.0.and.isol3.gt.0.and.isol4.gt.0) 
     .                                                         nphctt=10
C
       if(nphctt.eq.0) mulptt=1      ! to avoid unnecessary subdivisions
C
       GO TO 11
      ENDIF
C
C**** 6 NODED ELEMENT
C
      IF(NNODLT.EQ.6) THEN
       CALL RUNENDT('IDENPIT:ERROR DETECTED WITH NDIME=2')
       GO TO 11
      ENDIF
C
C**** 8 NODED ELEMENT
C
      IF(NNODLT.EQ.8) THEN
       CALL RUNENDT('IDENPIT:ERROR DETECTED WITH NDIME=2')
       GO TO 11
      ENDIF
C
C**** 9 NODED ELEMENT
C
      IF(NNODLT.EQ.9) THEN
       CALL RUNENDT('IDENPIT:ERROR DETECTED WITH NDIME=2')
       GO TO 11
      ENDIF
C
    3 CONTINUE
      CALL RUNENDT('IDENPIT:ERROR DETECTED WITH NDIME=3')
      GO TO 11
C
   11 CONTINUE
C
      IF(MULPDT.EQ.1) GO TO 22
      IF(NISDTT.EQ.2) GO TO 22
C
C**** IDENTIFY INTEGRATION ELEMENT FOR THE INTERNAL ENERGY FUNCTION IN
C     TIME t+dt FOR THE CASE NISDTT=1 (BIPHASE ELEMENT)
C
      NPHCDT=0
C
C**** ADDITIONAL CONTROL
C
      IF(NRULET.GT.5)
     . CALL RUNENDT('ERROR:ISOTH. PHASE-CHANGE WITH IMPROVED INT. RULE')
C
      GOTO(10,20,30), NDIMET
C
   10 CONTINUE
C
C**** 2 NODED ELEMENT
C
      IF(NNODLT.EQ.2) THEN
       NPHCDT=1
       GO TO 22
      ENDIF
C
C**** 3 NODED ELEMENT
C
      IF(NNODLT.EQ.3) THEN
       CALL RUNENDT('IDENPIT:ERROR DETECTED WITH NDIME=1')
c
c      not implemented yet
c
c      IF(ELDIST(1,1).GT.ELDIST(1,3).AND.ELDIST(1,2).GT.ELDIST(1,3))
c    .  NPHCDT=7
c      IF(ELDIST(1,1).LT.ELDIST(1,3).AND.ELDIST(1,2).LT.ELDIST(1,3))
c    .  NPHCDT=7
c      NPHCDT=6
       GO TO 22
      ENDIF
C
   20 CONTINUE
C
C**** 3 NODED ELEMENT
C
      IF(NNODLT.EQ.3) THEN
       ELDI1=ELDIST(1,1)
       ELDI2=ELDIST(1,2)
       ELDI3=ELDIST(1,3)
C
       CALL COIN02T(ELDI2,ELDI1,TEMDLT,RCOOR,ISOL1)
       CALL COIN02T(ELDI3,ELDI1,TEMDLT,RCOOR,ISOL2)
       CALL COIN02T(ELDI3,ELDI2,TEMDLT,RCOOR,ISOL3)
C
       if(isol1.gt.0.and.isol2.eq.0.and.isol3.gt.0) nphcdt=7
       if(isol1.gt.0.and.isol2.gt.0.and.isol3.eq.0) nphcdt=8
       if(isol1.eq.0.and.isol2.gt.0.and.isol3.gt.0) nphcdt=9
C
       if(nphcdt.eq.0) mulpdt=1      ! to avoid unnecessary subdivisions
C
       GO TO 22
      ENDIF
C
C**** 4 NODED ELEMENT
C
      IF(NNODLT.EQ.4) THEN
       ELDI1=ELDIST(1,1)
       ELDI2=ELDIST(1,2)
       ELDI3=ELDIST(1,3)
       ELDI4=ELDIST(1,4)
C
       CALL COIN01T(ELDI2,ELDI1,TEMDLT,RCOOR,ISOL1)
       CALL COIN01T(ELDI3,ELDI2,TEMDLT,RCOOR,ISOL2)
       CALL COIN01T(ELDI4,ELDI3,TEMDLT,RCOOR,ISOL3)
       CALL COIN01T(ELDI1,ELDI4,TEMDLT,RCOOR,ISOL4)
C
       if(isol1.gt.0.and.isol2.eq.0.and.isol3.eq.0.and.isol4.gt.0)
     .                                                         nphcdt=1
       if(isol1.gt.0.and.isol2.gt.0.and.isol3.eq.0.and.isol4.eq.0)
     .                                                         nphcdt=2
       if(isol1.eq.0.and.isol2.gt.0.and.isol3.gt.0.and.isol4.eq.0)
     .                                                         nphcdt=3
       if(isol1.eq.0.and.isol2.eq.0.and.isol3.gt.0.and.isol4.gt.0)
     .                                                         nphcdt=4
       if(isol1.gt.0.and.isol2.eq.0.and.isol3.gt.0.and.isol4.eq.0)
     .                                                         nphcdt=5
       if(isol1.eq.0.and.isol2.gt.0.and.isol3.eq.0.and.isol4.gt.0)
     .                                                         nphcdt=6
       if(isol1.gt.0.and.isol2.gt.0.and.isol3.gt.0.and.isol4.gt.0)
     .                                                         nphcdt=10
C
       if(nphcdt.eq.0) mulpdt=1      ! to avoid unnecessary subdivisions
C
       GO TO 22
      ENDIF
C
C**** 6 NODED ELEMENT
C
      IF(NNODLT.EQ.6) THEN
       CALL RUNENDT('IDENPIT:ERROR DETECTED WITH NDIME=2')
       GO TO 22
      ENDIF
C
C**** 8 NODED ELEMENT
C
      IF(NNODLT.EQ.8) THEN
       CALL RUNENDT('IDENPIT:ERROR DETECTED WITH NDIME=2')
       GO TO 22
      ENDIF
C
C**** 9 NODED ELEMENT
C
      IF(NNODLT.EQ.9) THEN
       CALL RUNENDT('IDENPIT:ERROR DETECTED WITH NDIME=2')
       GO TO 22
      ENDIF
C
   30 CONTINUE
      CALL RUNENDT('IDENPIT:ERROR DETECTED WITH NDIME=3')
      GO TO 22
C
   22 CONTINUE
C
      RETURN
      END
