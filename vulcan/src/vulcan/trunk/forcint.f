      SUBROUTINE FORCINT(DISTOT,ELDATT,ELPRET,ELVART,ELMATT,ELOADT,
     .                   HEADST,LNODST,MATNOT,PROELT,PROPST,WORK1T,
     .                   ELOATT,VELCTT,VEL1CT,DISPLT,PREAST,TGAPST,
     .                   COORDT,TEMPIT,ADVELT,BOUCHT,FPCHAT,LACTIT,
     .                   HTLODT)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE RESISTING FORCES
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'          ! thermal-mechanical
      INCLUDE 'nued_om.f'          ! thermal-microstructural
      INCLUDE 'nuef_om.f'          ! thermal-flow
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      COMMON/JACOBSTA/IERORT,KERORT
      COMMON/SMOGAUT/FACTAT
C
      DIMENSION MATNOT(NELEMT),        LNODST(NNODET,NELEMT),
     .          PROELT(NPRELT,NGRUPT), PROPST(NPROPT,NMATST),
     .          ELDATT(NDATAT),        ELPRET(NPREVT),
     .          ELVART(NSTATT),        ELMATT(NMATXT),
     .          WORK1T(*)
      DIMENSION DISTOT(NTOTVT,3),      ELOADT(*),
     .          ELOATT(*),             HEADST(*),
     .          VELCTT(*),             VEL1CT(*)
      DIMENSION DISPLT(NTOTVM)
      DIMENSION PREAST(*),             TGAPST(*)
      DIMENSION COORDT(NDIMET,NPOINT)
      DIMENSION TEMPIT(NPOINT),        ADVELT(NTOTVT*NDIMET),
     .          BOUCHT(NPOINT),        FPCHAT(NFPCH,NPOINT),
     .          LACTIT(NELEMT),        HTLODT(NHLODT,NSUBFT,NFUNCT)
C
      KUNLDT=0
      IF(LARGET.NE.0) KERORT=0
C
C**** LOOP OVER ELEMENTS
C
      DO 1000 IELEMT=1,NELEMT
      IF(NACTIT.EQ.1) THEN
       IF(LACTIT(IELEMT).EQ.0) GO TO 1000     ! skip non active elements
      ENDIF
      LGRUPT=MATNOT(IELEMT)
      LMATST=INT(PROELT(1,LGRUPT))
      NNODLT=INT(PROELT(2,LGRUPT))
      LTYPET=INT(PROELT(5,LGRUPT))
C
      NNODXT=NNODLT
      IF(LTYPET.EQ.104) THEN
       NOCOLT=INT(PROELT(11,LGRUPT))
       IF(NOCOLT.EQ.1) THEN                   ! non-coincident mesh
        NNODNT=INT(PROELT(10,LGRUPT))
        NNODXT=NNODNT
       ENDIF
      ENDIF
C
      IF(LARGET.NE.0) IERORT=0
C
C**** READ ELEM. DATA FROM DATA BASE
C
      IF(KPROBT.EQ.3) THEN
       CALL DATBAST(ELMATT(IMATXT(2)),    7,    2) ! ESTIF
      ELSE
       IF(NMEMO1.EQ.0.OR.NMEMO2.EQ.0)
     .  CALL DATBAST(ELDATT,    1,    2)   ! geometrical data
C
       IF(NMEMO.EQ.1)
     .  CALL DATBAST(ELPRET,    2,    2)   ! current values
C
       IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1.OR.NMEMO5.EQ.0)
     .  CALL DATBAST(ELVART,    5,    2)   ! last converged
C
       IF(NMEMO5.EQ.0) THEN
C
C**** SCALAR GATHER OPERATIONS ( DISTOT ---> ELDIST(ELVART) )
C
        CALL GATHER(DISTOT,NDOFCT,NPOINT,ELVART,         NDOFCT,NNODXT,
     .              LNODST(1,IELEMT))
       ELSE
C
C**** SCALAR GATHER OPERATIONS ( DISTOT ---> ELDIST )
C
        CALL GATHER(DISTOT,NDOFCT,NPOINT,WORK1T(IFORCT(44)),NDOFCT,
     .              NNODXT,
     .              LNODST(1,IELEMT))
       ENDIF                   ! nmemo5.eq.0
      END IF                   ! kprobt.eq.3
C
C**** SCALAR GATHER OPERATIONS ( VELCTT ---> WORK1T(IFORCT(12)) )
C
      IF(KDYNAT.EQ.1) THEN
       CALL GATHER(VELCTT,NDOFCT,NPOINT,WORK1T(IFORCT(12)),NDOFCT,
     .             NNODXT,
     .             LNODST(1,IELEMT))    ! current temp. rate
      ENDIF
C
C**** SCALAR GATHER OPERATIONS ( VEL1CT ---> WORK1T(IFORCT(33)) )
C
      IF(KDYNAT.EQ.1) THEN
       IF(IMICR.EQ.1) THEN
        CALL GATHER(VEL1CT,NDOFCT,NPOINT,WORK1T(IFORCT(33)),NDOFCT,
     .              NNODLT,
     .              LNODST(1,IELEMT))   ! predicted temp. rate
       ENDIF
      ENDIF
C
C**** SCALAR GATHER OPERATIONS ( TEMPIT ---> WORK1T(IFORCT(45)) )
C
      CALL GATHER(TEMPIT,NDOFCT,NPOINT,WORK1T(IFORCT(45)),NDOFCT,NNODLT,
     .            LNODST(1,IELEMT))     ! initial temperature
C
C**** SCALAR GATHER OPERATIONS ( BOUCHT ---> WORK1T(IFORCT(46)) )
C
      IF(NMEMO10.EQ.1)
     . CALL GATHER(BOUCHT,NDOFCT,NPOINT,WORK1T(IFORCT(46)),NDOFCT,
     .             NNODLT,
     .             LNODST(1,IELEMT))    ! boundary changes due to rho
C
      IF(ITERME.GT.0) THEN              ! bidirectional coupled
       IAUXX=0
       IF(ITERMG.GT.0) THEN             ! gap dependency
        IF(NMEMO11.EQ.0) THEN
C
C**** SCALAR GATHER OPERATIONS ( PREAST ---> WORK1T(IFORCT(14)) )
C
         CALL GATHER(PREAST,NDOFCT,NPOINT,WORK1T(IFORCT(14)),NDOFCT,
     .               NNODLT,
     .               LNODST(1,IELEMT))
C
C**** SCALAR GATHER OPERATIONS ( TGAPST ---> WORK1T(IFORCT(15)) )
C
         CALL GATHER(TGAPST,NDOFCT,NPOINT,WORK1T(IFORCT(15)),NDOFCT,
     .               NNODLT,
     .               LNODST(1,IELEMT))
        ELSE
C
C**** SCALAR GATHER OPERATIONS ( DISPLT ---> WORK1T(IFORCT(13)) )
C
         IAUXX=1
         CALL GATHER(DISPLT,NDOFCM,NPOINT,WORK1T(IFORCT(13)),NDOFCM,
     .               NNODLT,
     .               LNODST(1,IELEMT))
        ENDIF                       ! nmemo11.eq.0
       ENDIF                        ! itermg.gt.0
C
       IF((ITERMD.EQ.1.OR.ITERMF.EQ.1).AND.
     .                             IAUXX.EQ.0) THEN     ! deformed shape
C
C**** SCALAR GATHER OPERATIONS ( DISPLT ---> WORK1T(IFORCT(13)) )
C
        CALL GATHER(DISPLT,NDOFCM,NPOINT,WORK1T(IFORCT(13)),NDOFCM,
     .              NNODLT,
     .              LNODST(1,IELEMT))
       ENDIF                        ! itermd.eq.1. ... .and.iauxx.eq.0
      ENDIF                         ! iterme.gt.0
C
C**** SCALAR GATHER OPERATIONS ( ADVELT ---> WORK1T(IFORCT(34)) )
C
      IF(ICONVT.EQ.1)
     . CALL GATHER(ADVELT,NDIMET,NTOTVT,WORK1T(IFORCT(34)),NDIMET,
     .             NNODLT,
     .             LNODST(1,IELEMT))
C
C**** GATHER NODAL COORDINATES
C
      IF(NMEMO1.EQ.1) THEN    ! coordinates in a global array
       CALL GATHER(COORDT,NDIMET,NPOINT,WORK1T(IFORCT(16)),NDIMET,
     .             NNODLT,LNODST(1,IELEMT))
      ENDIF
C
C**** GATHER PHASE-CHANGE VECTOR
C
      IF(NFILL.EQ.1.OR.ITERMEF.GT.0)
     . CALL GATHER(FPCHAT,NFPCH,NPOINT,WORK1T(IFORCT(48)),NFPCH,
     .             NNODLT,LNODST(1,IELEMT))
C
C**** SCALAR GATHER OPERATIONS ( DISTOT ---> WORK1T(IFORCT(55)) )
C
      IF(ITERMEF.GT.0) THEN
       IF(KTEM1F.EQ.3) THEN
        NCOMPX=1+KDYNAT+NMEMO8+1
        CALL GATHER(DISTOT(1,NCOMPX),NDOFCT,NPOINT,WORK1T(IFORCT(55)),
     .              NDOFCT,NNODLT,
     .              LNODST(1,IELEMT)) ! coupling algorithmic temperature
       ENDIF
      ENDIF
C
C**** CALL ELEMENT PROCESSOR TO PERFORM REQUIRED OPERATIONS
C
      CALL ELMLIBT(LNODST(1,IELEMT),PROELT(1,LGRUPT),PROPST(1,LMATST),
     .             ELDATT,ELPRET,ELVART,ELMATT,HTLODT,WORK1T,     4)
C
      IF(LARGET.NE.0) THEN
       IF(IERORT.NE.0) GO TO 1000
      ENDIF
C
C**** WRITE CURRENT STATE VARIABLES TO DATA BASE
C 
      IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1.OR.NMEMO5.EQ.0)
     . CALL DATBAST(ELVART,    3,    1)
C
C**** SCALAR SCATTER OPERATION ( WORK1T --> ELOADT )
C
      CALL SCATER(WORK1T(IFORCT(1)), NDOFNT,NNODLT,ELOADT,NDOFCT,NPOINT,
     .            LNODST(1,IELEMT))
C
C**** SCALAR SCATTER OPERATION ( WORK1T --> ELOATT )
C
      CALL SCATER(WORK1T(IFORCT(11)),NDOFNT,NNODLT,ELOATT,NDOFCT,NPOINT,
     .            LNODST(1,IELEMT))
C
C**** SCALAR SCATTER OPERATION ( WORK1T --> BOUCHT )
C
      IF(NMEMO10.EQ.1)
     . CALL SCATERA(WORK1T(IFORCT(46)),NDOFCT,NNODLT,BOUCHT,NDOFCT,
     .              NPOINT,LNODST(1,IELEMT))
C
C**** SCALAR SCATTER OPERATION ( WORK1T --> FPCHAT )
C
      IF(ICONVT.EQ.1.AND.KDYNAT.EQ.0) THEN
       IF(LTYPET.EQ.5)
     .  CALL SCATERA(WORK1T(IFORCT(48)),NFPCH,NNODLT,FPCHAT,NFPCH,
     .               NPOINT,LNODST(1,IELEMT))
      ENDIF
C
 1000 CONTINUE
C
      IF(LARGET.NE.0) THEN
       IF(KERORT.NE.0)
     .  CALL RUNENDT('ERROR IN JACOBIAN MATRIX           ')
      ENDIF
C
C**** DEALS WITH ALTERNATIVE OPTION FOR MICRO ADVECTIVE PROBLEMS
C
      IF(ICONVT.EQ.1.AND.NMEMO3.EQ.1) THEN
       IF(IEGFPC.EQ.1) THEN     ! grad(f_pc) directly computed from f_pc
C
C**** COMPUTES L*f_pc & ITS RATE AT NODES
C
        FACTAT=1.D+00
C
        DO IPOINT=1,NPOINT
         WORK1T(IFORCT(51)+IPOINT-1)=0.D+00             ! accpn
         WORK1T(IFORCT(52)+IPOINT-1)=0.D+00             ! smstp
         WORK1T(IFORCT(52)+NPOINT+IPOINT-1)=0.D+00
        ENDDO
C
C**** LOOP OVER THE ELEMENTS
C
        DO 2000 IELEMT=1,NELEMT
        IF(NACTIT.EQ.1) THEN
         IF(LACTIT(IELEMT).EQ.0) GO TO 2000   ! skip non active elements
        ENDIF
        LGRUPT=MATNOT(IELEMT)
        LMATST=INT(PROELT(1,LGRUPT))
        LTYPET=INT(PROELT(5,LGRUPT))
        NNODLT=INT(PROELT(2,LGRUPT))
        IF(LTYPET.EQ.101.OR.LTYPET.EQ.104) GO TO 2000 ! SKIP CONTACT EL.
C
C**** READ ELDAT & ELVAR FROM DATA BASE
C
        IF(NMEMO1.EQ.0.OR.NMEMO2.EQ.0)
     .   CALL DATBAST(ELDATT,    1,    2)
C
        IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1.OR.NMEMO5.EQ.0)
     .   CALL DATBAST(ELVART,    3,    2)    ! current
C
C**** GATHER NODAL COORDINATES
C
        IF(NMEMO1.EQ.1) THEN    ! coordinates in a global array
         CALL GATHER(COORDT,NDIMET,NPOINT,WORK1T(IFORCT(16)),NDIMET,
     .               NNODLT,LNODST(1,IELEMT))
        ENDIF
C
        IF(ITERME.GT.0) THEN          ! bidirectional coupled
         IF(ITERMD.EQ.1) THEN         ! deformed shape
C
C**** SCALAR GATHER OPERATIONS ( DISPLT ---> WORK1T(IFORCT(13)) )
C
          CALL GATHER(DISPLT,NDOFCM,NPOINT,WORK1T(IFORCT(13)),NDOFCM,
     .                NNODLT,
     .                LNODST(1,IELEMT))
         ENDIF                        ! itermd.eq.1
        ENDIF                         ! iterme.gt.0
C
C**** COMPUTES ACCPNT
C
        DO INODET=1,NNODST
         LPOINT=LNODST(INODET,IELEMT)
         IF(LPOINT.NE.0) WORK1T(IFORCT(51)+LPOINT-1)=
     .                                WORK1T(IFORCT(51)+LPOINT-1)+FACTAT
        ENDDO
C
C**** CALL THE ELEMENT PROCESSOR TO EVALUATE ELEMENT CONTRIBUTION TO
C     GLOBAL VECTOR OF L*f_pc & ITS RATE
C
        CALL ELMLIBT(LNODST(1,IELEMT),PROELT(1,LGRUPT),PROPST(1,LMATST),
     .               ELDATT,ELPRET,ELVART,ELMATT,DUMMYT,WORK1T,    19)
C
C**** INCLUDE IT IN A GLOBAL VECTOR
C
        CALL OUTASS(LNODST,WORK1T(IFORCT(52)),WORK1T(IFORCT(54)),FACTAT,
     .              NNODET,IELEMT,NELEMT,NPOINT,
     .                   2)
C
 2000   CONTINUE
C
C**** PERFORMS AVERAGE MEAN FOR L*f_pc & ITS RATE (KDYNAT=1 is assumed)
C
        IFPCH=2*NNUPT+NNUPO+NFILL+IGALFA+1
        DO IPOINT=1,NPOINT
         IF(WORK1T(IFORCT(51)+IPOINT-1).NE.0.0D0) THEN
          FPCHAT(IFPCH  ,IPOINT)=WORK1T(IFORCT(52)+IPOINT-1)/
     .                                       WORK1T(IFORCT(51)+IPOINT-1)
          FPCHAT(IFPCH+1,IPOINT)=
     .                              (WORK1T(IFORCT(52)+NPOINT+IPOINT-1)/
     .                               WORK1T(IFORCT(51)+IPOINT-1))/DTIMET
         ENDIF
        ENDDO
C
C**** INTERPOLATE FOR MID-SIDE NODES
C
        DO 3000 IELEMT=1,NELEMT
        IF(NACTIT.EQ.1) THEN
         IF(LACTIT(IELEMT).EQ.0) GO TO 3000   ! skip non active elements
        ENDIF
        LGRUPT=MATNOT(IELEMT)
        LMATST=INT(PROELT(1,LGRUPT))
        LTYPET=INT(PROELT(5,LGRUPT))
        NNODLT=INT(PROELT(2,LGRUPT))
        IF(LTYPET.EQ.101.OR.LTYPET.EQ.104) GO TO 3000 ! SKIP CONTACT EL.
C
C****  CALL THE ELEMENT PROCESSOR TO NOTHING (computes NQUTRT)
C
        CALL ELMLIBT(LNODST(1,IELEMT),PROELT(1,LGRUPT),PROPST(1,LMATST),
     .               ELDATT,ELPRET,ELVART,ELMATT,DUMMYT,WORK1T,    18)
C
        CALL SMONOD(NDIMET,NNODLT,NNODST,NQUTRT)
C
        IF(NNODLT.EQ.NNODST) GO TO 3000
C
C**** L*f_pc
C
        DO INODET=1,NNODST
         IPOINT=LNODST(INODET,IELEMT)
         WORK1T(IFORCT(53)+INODET-1)=WORK1T(IFORCT(52)+IPOINT-1)
        END DO
        CALL SMOMID(WORK1T(IFORCT(53)),NDIMET,NNODLT,NQUTRT)
        DO INODET=NNODST+1,NNODLT
         IPOINT=LNODST(INODET,IELEMT)
         FPCHAT(IFPCH,IPOINT)=WORK1T(IFORCT(53)+INODET-1)
        END DO
C
C**** RATE OF L*f_pc
C
        IF(KDYNAT.EQ.1) THEN
         DO INODET=1,NNODST
          IPOINT=LNODST(INODET,IELEMT)
          WORK1T(IFORCT(53)+INODET-1)=WORK1T(IFORCT(52)+NPOINT+IPOINT-1)
         END DO
         CALL SMOMID(WORK1T(IFORCT(53)),NDIMET,NNODLT,NQUTRT)
         DO INODET=NNODST+1,NNODLT
          IPOINT=LNODST(INODET,IELEMT)
          FPCHAT(IFPCH+1,IPOINT)=WORK1T(IFORCT(53)+INODET-1)
         END DO
        ENDIF
C
 3000   CONTINUE
C
C**** LOOP OVER THE ELEMENTS (computes residual cont. of ad-pc term)
C
        DO 4000 IELEMT=1,NELEMT
        IF(NACTIT.EQ.1) THEN
         IF(LACTIT(IELEMT).EQ.0) GO TO 4000   ! skip non active elements
        ENDIF
        LGRUPT=MATNOT(IELEMT)
        LMATST=INT(PROELT(1,LGRUPT))
        NNODLT=INT(PROELT(2,LGRUPT))
        LTYPET=INT(PROELT(5,LGRUPT))
C
C**** READ ELEM. DATA FROM DATA BASE
C
        IF(KPROBT.EQ.3) THEN
         CALL DATBAST(ELMATT(IMATXT(2)),    7,    2) ! ESTIF
        ELSE
         IF(NMEMO1.EQ.0.OR.NMEMO2.EQ.0)
     .    CALL DATBAST(ELDATT,    1,    2)   ! geometrical data
C
         IF(NMEMO.EQ.1)
     .    CALL DATBAST(ELPRET,    2,    2)   ! current values
C
         IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1.OR.NMEMO5.EQ.0)
     .    CALL DATBAST(ELVART,    3,    2)   ! current values
C
         IF(NMEMO5.EQ.0) THEN
C
C**** SCALAR GATHER OPERATIONS ( DISTOT ---> ELDIST(ELVART) )
C
          CALL GATHER(DISTOT,NDOFCT,NPOINT,ELVART,            NDOFCT,
     .                NNODLT,
     .                LNODST(1,IELEMT))
         ELSE
C
C**** SCALAR GATHER OPERATIONS ( DISTOT ---> ELDIST )
C
          CALL GATHER(DISTOT,NDOFCT,NPOINT,WORK1T(IFORCT(44)),NDOFCT,
     .                NNODLT,
     .                LNODST(1,IELEMT))
         ENDIF                   ! nmemo5.eq.0
        END IF                   ! kprobt.eq.3
C
C**** SCALAR GATHER OPERATIONS ( VELCTT ---> WORK1T(IFORCT(12)) )
C
        IF(KDYNAT.EQ.1) THEN
         CALL GATHER(VELCTT,NDOFCT,NPOINT,WORK1T(IFORCT(12)),NDOFCT,
     .               NNODLT,
     .               LNODST(1,IELEMT))    ! current temp. rate
        ENDIF
C
C**** SCALAR GATHER OPERATIONS ( VEL1CT ---> WORK1T(IFORCT(33)) )
C
        IF(KDYNAT.EQ.1) THEN
         IF(IMICR.EQ.1) THEN
          CALL GATHER(VEL1CT,NDOFCT,NPOINT,WORK1T(IFORCT(33)),NDOFCT,
     .                NNODLT,
     .                LNODST(1,IELEMT))   ! predicted temp. rate
         ENDIF
        ENDIF
C
C**** SCALAR GATHER OPERATIONS ( TEMPIT ---> WORK1T(IFORCT(45)) )
C
        CALL GATHER(TEMPIT,NDOFCT,NPOINT,WORK1T(IFORCT(45)),NDOFCT,
     .              NNODLT,
     .              LNODST(1,IELEMT))     ! initial temperature
C
        IF(ITERME.GT.0) THEN              ! bidirectional coupled
         IAUXX=0
         IF(ITERMG.GT.0) THEN             ! gap dependency
          IF(NMEMO11.EQ.0) THEN
C
C**** SCALAR GATHER OPERATIONS ( PREAST ---> WORK1T(IFORCT(14)) )
C
           CALL GATHER(PREAST,NDOFCT,NPOINT,WORK1T(IFORCT(14)),NDOFCT,
     .                 NNODLT,
     .                 LNODST(1,IELEMT))
C
C**** SCALAR GATHER OPERATIONS ( TGAPST ---> WORK1T(IFORCT(15)) )
C
           CALL GATHER(TGAPST,NDOFCT,NPOINT,WORK1T(IFORCT(15)),NDOFCT,
     .                 NNODLT,
     .               LNODST(1,IELEMT))
          ELSE
C
C**** SCALAR GATHER OPERATIONS ( DISPLT ---> WORK1T(IFORCT(13)) )
C
           IAUXX=1
           CALL GATHER(DISPLT,NDOFCM,NPOINT,WORK1T(IFORCT(13)),NDOFCM,
     .                 NNODLT,
     .                 LNODST(1,IELEMT))
          ENDIF                       ! nmemo11.eq.0
         ENDIF                        ! itermg.gt.0
C
         IF((ITERMD.EQ.1.OR.ITERMF.EQ.1).AND.
     .                             IAUXX.EQ.0) THEN     ! deformed shape
C
C**** SCALAR GATHER OPERATIONS ( DISPLT ---> WORK1T(IFORCT(13)) )
C
          CALL GATHER(DISPLT,NDOFCM,NPOINT,WORK1T(IFORCT(13)),NDOFCM,
     .                NNODLT,
     .                LNODST(1,IELEMT))
         ENDIF                        ! itermd.eq.1. ... .and.iauxx.eq.0
        ENDIF                         ! iterme.gt.0
C
C**** SCALAR GATHER OPERATIONS ( ADVELT ---> WORK1T(IFORCT(34)) )
C
        IF(ICONVT.EQ.1)
     .   CALL GATHER(ADVELT,NDIMET,NTOTVT,WORK1T(IFORCT(34)),NDIMET,
     .               NNODLT,
     .               LNODST(1,IELEMT))
C
C**** GATHER NODAL COORDINATES
C
        IF(NMEMO1.EQ.1) THEN    ! coordinates in a global array
         CALL GATHER(COORDT,NDIMET,NPOINT,WORK1T(IFORCT(16)),NDIMET,
     .               NNODLT,LNODST(1,IELEMT))
        ENDIF
C
C**** GATHER PHASE-CHANGE VECTOR
C
        IF(NFILL.EQ.1.OR.ITERMEF.GT.0)
     .   CALL GATHER(FPCHAT,NFPCH,NPOINT,WORK1T(IFORCT(48)),NFPCH,
     .               NNODLT,LNODST(1,IELEMT))
C
C**** CALL ELEMENT PROCESSOR TO PERFORM REQUIRED OPERATIONS
C
        CALL ELMLIBT(LNODST(1,IELEMT),PROELT(1,LGRUPT),PROPST(1,LMATST),
     .               ELDATT,ELPRET,ELVART,ELMATT,DUMMYT,WORK1T,     5)
C
C**** SCALAR SCATTER OPERATION ( WORK1T --> ELOADT )
C
        CALL SCATER(WORK1T(IFORCT(1)), NDOFNT,NNODLT,ELOADT,NDOFCT,
     .              NPOINT,
     .              LNODST(1,IELEMT))
C
 4000   CONTINUE
C
       ENDIF       ! iegfpc.eq.1
      ENDIF        ! iconvt.eq.1.and.imicr.eq.1
C
      RETURN
      END
