      SUBROUTINE OUTSMOT(ELDATT,ELPRET,ELVART,ELMATT,LNODST,MATNOT,
     .                   PROELT,PROPST,COORDT,PREAST,FPCHAT,DISTOT,
     .                   VELCTT,DISPLT,LACTIT,
     .                   SMSTST,SMSTNT,ACCPNT,SMSTPT,SMSPPT,
     .                   WORK1T)
C***********************************************************************
C
C**** THIS ROUTINE PERFORMS A SMOOTHING OVER THE HEAT FLUXES, TEMPERAT.
C     GRADIENTS & INTERNAL VARIABLES
C
C     MOREOVER, THIS ROUTINE ASSIGNS MICROSCOPICAL PHASE-CHANGE
C     FUNCTIONS TO FPCHAT
C
C
C     Notes:
C
C     The call to OUTASS includes NNODET (instead of NNODST) because the
C     dimension of the arrays in this routine are fixed
C
C***********************************************************************
C
C     Index of variables:
C
C     NNUPT: total number of phase-changes (maximum)
C     NNUPC: number of macroscopic phase-changes (maximum)
C     NNUPM: number of microscopic phase-changes (maximum)
C     NNUPO: other microst. variables to transfer to mechanical problem
C            by means of FPHCAT
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
      INCLUDE 'nued_om.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inpo_omt.f'
C
      COMMON/SMOGAUT/FACTAT
      COMMON/POROSI1/KPOROA
C
      DIMENSION ACCPNT(NPOINT),        SMSTST(NSTR1T,NPOINT),
     .          SMSTNT(NSTR1T,NPOINT), SMSTPT(NNUINT,NPOINT),
     .          SMSPPT(NPOROT,NPOINT)
      DIMENSION MATNOT(NELEMT),        LNODST(NNODET,NELEMT),
     .          PROELT(NPRELT,NGRUPT), PROPST(NPROPT,NMATST),
     .          ELDATT(NDATAT),        ELPRET(NPREVT),
     .          ELVART(NSTATT),        ELMATT(NMATXT)
      DIMENSION COORDT(NDIMET,NPOINT), PREAST(NPOROT,NPOINT),
     .          FPCHAT(NFPCH,NPOINT),  DISTOT(*),
     .          VELCTT(*),             DISPLT(NTOTVM),
     .          LACTIT(NELEMT)
      DIMENSION WORK1T(*)
C
      FACTAT=1.D+00
      KPOROA=0
C
      DO IPOINT=1,NPOINT
       ACCPNT(IPOINT)=0.D+00
       IF(NMEMO4.EQ.1) THEN
        DO ISTRET=1,NSTR1T
         SMSTST(ISTRET,IPOINT)=0.D+00
         SMSTNT(ISTRET,IPOINT)=0.D+00
        ENDDO
       ENDIF
       IF(NNUINT.GT.0) THEN
        DO INUINT=1,NNUINT
         SMSTPT(INUINT,IPOINT)=0.D+00
        ENDDO
       ENDIF
       IF(NPOROT.GT.0) THEN       ! ICALPO not considered; see outputt.f
        DO IPOROT=1,NPOROT
         SMSPPT(IPOROT,IPOINT)=0.D+00
        ENDDO
       ENDIF
      ENDDO
C
C**** LOOP OVER THE ELEMENTS
C
      DO 100 IELEMT=1,NELEMT
       IF(NACTIT.EQ.1) THEN
        IF(LACTIT(IELEMT).EQ.0) GO TO 100     ! skip non active elements
       ENDIF
       LGRUPT=MATNOT(IELEMT)
       LMATST=INT(PROELT(1,LGRUPT))
       LTYPET=INT(PROELT(5,LGRUPT))
       NNODLT=INT(PROELT(2,LGRUPT))
       IF(LTYPET.EQ.101.OR.LTYPET.EQ.104) GO TO 100  !SKIP CONTACT ELEM.
C
C**** READ ELDAT & ELVAR FROM DATA BASE
C
       IF(NMEMO1.EQ.0.OR.NMEMO2.EQ.0)
     .  CALL DATBAST(ELDATT,    1,    2)
C
       IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1.OR.NMEMO5.EQ.0)
     .  CALL DATBAST(ELVART,    3,    2)    ! current
C
C**** GATHER NODAL COORDINATES
C
       IF(NMEMO1.EQ.1) THEN    ! coordinates in a global array
        CALL GATHER(COORDT,NDIMET,NPOINT,WORK1T(IGSMOT(19)),NDIMET,
     .              NNODLT,LNODST(1,IELEMT))
       ENDIF
C
       IF(NMEMO4.EQ.1) THEN
        IF(ITERME.GT.0) THEN          ! bidirectional coupled
         IF(ITERMD.EQ.1) THEN         ! deformed shape
C
C**** SCALAR GATHER OPERATIONS ( DISPLT ---> WORK1T(IGSMOT(30)) )
C
          CALL GATHER(DISPLT,NDOFCM,NPOINT,WORK1T(IGSMOT(30)),NDOFCM,
     .                NNODLT,
     .                LNODST(1,IELEMT))
         ENDIF                        ! itermd.eq.1
        ENDIF                         ! iterme.gt.0
C
C**** CALL THE ELEMENT PROCESSOR TO EVALUATE ELEMENT CONTRIBUTION TO
C     GLOBAL VECTOR OF HEAT FLUXES
C
        CALL ELMLIBT(LNODST(1,IELEMT),PROELT(1,LGRUPT),PROPST(1,LMATST),
     .               ELDATT,ELPRET,ELVART,ELMATT,DUMMYT,WORK1T,    13)
C
C**** COMPUTES ACCPNT
C
        DO INODET=1,NNODST
         LPOINT=LNODST(INODET,IELEMT)
         IF(LPOINT.NE.0) ACCPNT(LPOINT)=ACCPNT(LPOINT)+FACTAT
        ENDDO
C
C**** INCLUDE IT IN A GLOBAL VECTOR
C
        CALL OUTASS(LNODST,SMSTST,WORK1T(IGSMOT(5)),FACTAT,
     .              NNODET,IELEMT,NELEMT,NPOINT,
     .              NSTR1T)
C
C**** CALL THE ELEMENT PROCESSOR TO EVALUATE ELEMENT CONTRIBUTION TO
C     GLOBAL VECTOR OF TEMPERATURE GRADIENTS (not implemented)
C
c       CALL ELMLIBT(LNODST(1,IELEMT),PROELT(1,LGRUPT),PROPST(1,LMATST),
c    .               ELDATT,ELPRET,ELVART,ELMATT,DUMMYT,WORK1T,    14)
C
C**** INCLUDE IT IN A GLOBAL VECTOR
C
c       CALL OUTASS(LNODST,SMSTNT,WORK1T(IGSMOT(7)),FACTAT,
c    .              NNODET,IELEMT,NELEMT,NPOINT,
c    .              NSTR1T)
       ENDIF     ! nmemo4.eq.1
C
       IF(NNUINT.GT.0) THEN
        IF(NMEMO4.EQ.0) THEN
         IF(ITERME.GT.0) THEN          ! bidirectional coupled
          IF(ITERMD.EQ.1) THEN         ! deformed shape
C
C**** SCALAR GATHER OPERATIONS ( DISPLT ---> WORK1T(IGSMOT(30)) )
C
           CALL GATHER(DISPLT,NDOFCM,NPOINT,WORK1T(IGSMOT(30)),NDOFCM,
     .                 NNODLT,
     .                 LNODST(1,IELEMT))
          ENDIF                        ! itermd.eq.1
         ENDIF                         ! iterme.gt.0
        ENDIF                          ! nmemo4.eq.0
C
C**** CALL THE ELEMENT PROCESSOR TO EVALUATE ELEMENT CONTRIBUTION TO
C     GLOBAL VECTOR OF INTERNAL VARIABLES
C
        CALL ELMLIBT(LNODST(1,IELEMT),PROELT(1,LGRUPT),PROPST(1,LMATST),
     .               ELDATT,ELPRET,ELVART,ELMATT,DUMMYT,WORK1T,    15)
C
C**** COMPUTES ACCPN (IF NECESSARY)
C
        IF(NMEMO4.EQ.0) THEN
         DO INODET=1,NNODST
          LPOINT=LNODST(INODET,IELEMT)
          IF(LPOINT.NE.0) ACCPNT(LPOINT)=ACCPNT(LPOINT)+FACTAT
         ENDDO
        ENDIF
C
C**** INCLUDE IT IN A GLOBAL VECTOR
C
        CALL OUTASS(LNODST,SMSTPT,WORK1T(IGSMOT(10)),FACTAT,
     .              NNODET,IELEMT,NELEMT,NPOINT,
     .              NNUINT)
C
       ENDIF     ! nnuint.gt.0
C
C**** EVALUATE ELEM. CONTRIBUTION TO GLOBAL VECTOR OF POROSITY CRITERIA
C
       IF(NPOROT.GT.0.AND.ICALPO.EQ.1) THEN
        IF(IMICR.EQ.0) THEN
C
C**** GATHER NODAL POROSITY CRITERIA
C
         NGAULT=INT(PROELT(4,LGRUPT))
         NNODGT=NNODLT
         IF(NGAULT.GT.NNODLT) NNODGT=NGAULT
         CALL GATHER(PREAST,NPOROT,NPOINT,WORK1T(IGSMOT(25)),NPOROT,
     .               NNODGT,LNODST(1,IELEMT))
C
         IF(NMEMO5.EQ.0) THEN
C
C**** SCALAR GATHER OPERATIONS ( DISTOT ---> ELDIST(ELVART) )
C
          CALL GATHER(DISTOT,NDOFCT,NPOINT,ELVART,NDOFCT,NNODLT,
     .                LNODST(1,IELEMT))
         ELSE
C
C**** SCALAR GATHER OPERATIONS ( DISTOT ---> ELDIST )
C
          CALL GATHER(DISTOT,NDOFCT,NPOINT,WORK1T(IGSMOT(21)),NDOFCT,
     .                NNODLT,
     .                LNODST(1,IELEMT))
         ENDIF                   ! nmemo5.eq.0
C
C**** GATHER NODAL PHASE-CHANGE FUNCTION (MACROSCOPIC)
C
         CALL GATHER(FPCHAT,NFPCH,NPOINT,WORK1T(IGSMOT(26)),NFPCH,
     .               NNODLT,LNODST(1,IELEMT))
C
C**** SCALAR GATHER OPERATIONS ( VELCTT ---> VELCMT ) - TEMP. RATE
C
         IF(KDYNAT.EQ.1)
     .    CALL GATHER(VELCTT,NDOFCT,NPOINT,WORK1T(IGSMOT(27)),NDOFCT,
     .                NNODLT,
     .                LNODST(1,IELEMT))
C
         IF(NMEMO4.EQ.0.AND.NNUINT.EQ.0) THEN
          IF(ITERME.GT.0) THEN          ! bidirectional coupled
           IF(ITERMD.EQ.1) THEN         ! deformed shape
C
C**** SCALAR GATHER OPERATIONS ( DISPLT ---> WORK1T(IGSMOT(30)) )
C
            CALL GATHER(DISPLT,NDOFCM,NPOINT,WORK1T(IGSMOT(30)),NDOFCM,
     .                  NNODLT,
     .                  LNODST(1,IELEMT))
           ENDIF                        ! itermd.eq.1
          ENDIF                         ! iterme.gt.0
         ENDIF                          ! nmemo4.eq.0.and.nnuint.eq.0
C
C**** CALL THE ELEMENT PROCESSOR TO EVALUATE ELEMENT CONTRIBUTION TO
C     GLOBAL VECTOR OF POROSITY CRITERIA
C
         CALL ELMLIBT(LNODST(1,IELEMT),PROELT(1,LGRUPT),
     .                PROPST(1,LMATST),
     .                ELDATT,ELPRET,ELVART,ELMATT,DUMMYT,WORK1T,    17)
C
C**** ADDITIONAL CONTROL
C
         KPOROA=KPOROA+KPOROT
C
C**** COMPUTES ACCPN (IF NECESSARY)
C
         IF(NMEMO4.EQ.0.AND.NNUINT.EQ.0) THEN
          DO INODET=1,NNODST
           LPOINT=LNODST(INODET,IELEMT)
           IF(LPOINT.NE.0) ACCPNT(LPOINT)=ACCPNT(LPOINT)+FACTAT
          ENDDO
         ENDIF
C
C**** INCLUDE IT IN A GLOBAL VECTOR
C
         CALL OUTASS(LNODST,SMSPPT,WORK1T(IGSMOT(24)),FACTAT,
     .               NNODET,IELEMT,NELEMT,NPOINT,
     .               NPOROT)
C
        ENDIF    ! imicr.eq.0
       ENDIF     ! nporot.gt.0.and.icalpo.eq.1
C
  100 CONTINUE
C
C**** ASSIGNS MICROSCOPICAL PHASE-CHANGE FUNCTION TO FPCHAT
C
C     Notes:
C
C     IPLUAT(INUPM) (defined in pointes.f and passed through auxl_omt.f)
C     gives the position of the phase-change functions in
C     SMSTPT(IPLUAT(INUPM),IPOINT)
C
C     FPCHAT can be transferred to mechanical problem (coupled problems)
C     or be used in the porosity criteria (uncoupled & coupled problems)
C     or be transferred to flow problem (thermally coupled flows)
C
      IF(IMICR.EQ.1) THEN
       IF(NNUPM.GT.0) THEN
        DO INUPM=1,NNUPM
         INUAU=IPLUAT(INUPM)
         DO IPOINT=1,NPOINT
          FAUXX=FPCHAT(INUPM+NNUPC,IPOINT)              ! last converged
          IF(ACCPNT(IPOINT).NE.0.0D0)
     .     FPCHAT(INUPM+NNUPC,IPOINT)=SMSTPT(INUAU,IPOINT)/    ! current
     .                                                    ACCPNT(IPOINT)
C
          IF(FPCHAT(INUPM+NNUPC,IPOINT).GT.1.0D0)       ! control
     .                                  FPCHAT(INUPM+NNUPC,IPOINT)=1.0D0
          IF(FPCHAT(INUPM+NNUPC,IPOINT).LT.0.0D0)
     .                                  FPCHAT(INUPM+NNUPC,IPOINT)=0.0D0
C
          IF(IPRCOT.GT.0)
     .     FPCHAT(INUPM+NNUPC+NNUPT,IPOINT)=
     .     (FPCHAT(INUPM+NNUPC,IPOINT)-FAUXX)/DTIMET          ! dot f_pc
         ENDDO
        ENDDO
       ENDIF            ! nnupm.gt.0
C
C**** ASSIGNS OTHER VARIABLES TO FPCHAT
C
C     Notes:
C
C     IPLUOT(INUPO) (defined in pointes.f and passed through auxl_omt.f)
C     gives the position of the other variables in
C     SMSTPT(IPLUOT(INUPO),IPOINT)
C
C     FPCHAT can be transferred to mechanical problem (coupled problems)
C     or be used in the porosity criteria (uncoupled & coupled problems)
C     or be used in the thermal problem with advective effects when
C     solving the microstructural evolution equations (thermally
C     uncoupled & coupled problems)
C
       IF(NNUPO.GT.0) THEN
        DO INUPO=1,NNUPO
         INUAU=IPLUOT(INUPO)
         DO IPOINT=1,NPOINT
          IF(ACCPNT(IPOINT).NE.0.0D0)
     .     FPCHAT(INUPO+2*NNUPT,IPOINT)=SMSTPT(INUAU,IPOINT)/
     .                                                    ACCPNT(IPOINT)
         ENDDO
        ENDDO
       ENDIF            ! nnupo.gt.0
C
C**** LOOP OVER THE ELEMENTS
C
       DO 101 IELEMT=1,NELEMT
        IF(NACTIT.EQ.1) THEN
         IF(LACTIT(IELEMT).EQ.0) GO TO 101    ! skip non active elements
        ENDIF
        LGRUPT=MATNOT(IELEMT)
        LMATST=INT(PROELT(1,LGRUPT))
        LTYPET=INT(PROELT(5,LGRUPT))
        NNODLT=INT(PROELT(2,LGRUPT))
        IF(LTYPET.EQ.101.OR.LTYPET.EQ.104) GO TO 101 !SKIP CONTACT ELEM.
C
C**** READ ELDAT & ELVAR FROM DATA BASE
C
        IF(NMEMO1.EQ.0.OR.NMEMO2.EQ.0)
     .   CALL DATBAST(ELDATT,    1,    2)
C
        IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1)
     .   CALL DATBAST(ELVART,    3,    2)    ! current
C
C**** GATHER NODAL COORDINATES
C
        IF(NMEMO1.EQ.1) THEN    ! coordinates in a global array
         CALL GATHER(COORDT,NDIMET,NPOINT,WORK1T(IGSMOT(19)),NDIMET,
     .               NNODLT,LNODST(1,IELEMT))
        ENDIF
C
C**** EVALUATE ELEM. CONTRIBUTION TO GLOBAL VECTOR OF POROSITY CRITERIA
C
C     Note: the porosity criteria may contain microstructural
C           phase-change functions and other microstructural variables
C
        IF(NPOROT.GT.0.AND.ICALPO.EQ.1) THEN
C
C**** GATHER NODAL POROSITY CRITERIA
C
         NGAULT=INT(PROELT(4,LGRUPT))
         NNODGT=NNODLT
         IF(NGAULT.GT.NNODLT) NNODGT=NGAULT
         CALL GATHER(PREAST,NPOROT,NPOINT,WORK1T(IGSMOT(25)),NPOROT,
     .               NNODGT,LNODST(1,IELEMT))
C
         IF(NMEMO5.EQ.0) THEN
C
C**** SCALAR GATHER OPERATIONS ( DISTOT ---> ELDIST(ELVART) )
C
          CALL GATHER(DISTOT,NDOFCT,NPOINT,ELVART,NDOFCT,NNODLT,
     .                LNODST(1,IELEMT))
         ELSE
C
C**** SCALAR GATHER OPERATIONS ( DISTOT ---> ELDIST )
C
          CALL GATHER(DISTOT,NDOFCT,NPOINT,WORK1T(IGSMOT(21)),NDOFCT,
     .                NNODLT,
     .                LNODST(1,IELEMT))
         ENDIF                   ! nmemo5.eq.0
C
C**** GATHER NODAL PHASE-CHANGE FUNCTION (MACROSCOPIC, MICROSCOPIC &
C     OTHER MICROSTRUCTURE VARIABLES)
C
         CALL GATHER(FPCHAT,NFPCH,NPOINT,WORK1T(IGSMOT(26)),NFPCH,
     .               NNODLT,LNODST(1,IELEMT))
C
C**** SCALAR GATHER OPERATIONS ( VELCTT ---> VELCMT ) - TEMP. RATE
C
         IF(KDYNAT.EQ.1)
     .    CALL GATHER(VELCTT,NDOFCT,NPOINT,WORK1T(IGSMOT(27)),NDOFCT,
     .                NNODLT,
     .                LNODST(1,IELEMT))
C
         IF(NMEMO4.EQ.0.AND.NNUINT.EQ.0) THEN
          IF(ITERME.GT.0) THEN          ! bidirectional coupled
           IF(ITERMD.EQ.1) THEN         ! deformed shape
C
C**** SCALAR GATHER OPERATIONS ( DISPLT ---> WORK1T(IGSMOT(30)) )
C
            CALL GATHER(DISPLT,NDOFCM,NPOINT,WORK1T(IGSMOT(30)),NDOFCM,
     .                  NNODLT,
     .                  LNODST(1,IELEMT))
           ENDIF                        ! itermd.eq.1
          ENDIF                         ! iterme.gt.0
         ENDIF                          ! nmemo4.eq.0.and.nnuint.eq.0
C
C**** CALL THE ELEMENT PROCESSOR TO EVALUATE ELEMENT CONTRIBUTION TO
C     GLOBAL VECTOR OF POROSITY CRITERIA
C
         CALL ELMLIBT(LNODST(1,IELEMT),PROELT(1,LGRUPT),
     .                PROPST(1,LMATST),
     .                ELDATT,ELPRET,ELVART,ELMATT,DUMMYT,WORK1T,    17)
C
C**** ADDITIONAL CONTROL
C
         KPOROA=KPOROA+KPOROT
C
C**** COMPUTES ACCPN (IF NECESSARY)
C
         IF(NMEMO4.EQ.0.AND.NNUINT.EQ.0) THEN
          DO INODET=1,NNODST
           LPOINT=LNODST(INODET,IELEMT)
           IF(LPOINT.NE.0) ACCPNT(LPOINT)=ACCPNT(LPOINT)+FACTAT
          ENDDO
         ENDIF
C
C**** INCLUDE IT IN A GLOBAL VECTOR
C
         CALL OUTASS(LNODST,SMSPPT,WORK1T(IGSMOT(24)),FACTAT,
     .               NNODET,IELEMT,NELEMT,NPOINT,
     .               NPOROT)
C
        ENDIF     ! nporot.gt.0.and.icalpo.eq.1
C
  101  CONTINUE
      ENDIF             ! imicr.eq.1
C
C**** PERFORMS AVERAGE MEAN FOR HEAT FLUXES
C
      IF(NMEMO4.EQ.1) THEN
       DO IPOINT=1,NPOINT
        DO ISTRET=1,NSTR1T
         IF(ACCPNT(IPOINT).NE.0.0D0)
     .    SMSTST(ISTRET,IPOINT)=SMSTST(ISTRET,IPOINT)/ACCPNT(IPOINT)
        ENDDO
       ENDDO
      ENDIF
C
C**** PERFORMS AVERAGE MEAN FOR TEMPERATURE GRADIENTS (not implem. yet)
C

C
C**** PERFORMS AVERAGE MEAN FOR INTERNAL VARIABLES
C
      IF(NNUINT.GT.0) THEN
       DO IPOINT=1,NPOINT
        DO INUINT=1,NNUINT
         IF(ACCPNT(IPOINT).NE.0.0D0)
     .    SMSTPT(INUINT,IPOINT)=SMSTPT(INUINT,IPOINT)/ACCPNT(IPOINT)
        ENDDO
       ENDDO
      ENDIF
C
C**** PERFORMS AVERAGE MEAN FOR POROSITY CRITERIA & ASSIGNS TO PREAST
C
      IF(NPOROT.GT.0.AND.ICALPO.EQ.1) THEN
       DO IPOINT=1,NPOINT
        DO IPOROT=1,NPOROT
         IF(ACCPNT(IPOINT).NE.0.0D0) THEN
          SMSPPT(IPOROT,IPOINT)=SMSPPT(IPOROT,IPOINT)/ACCPNT(IPOINT)
          PREAST(IPOROT,IPOINT)=SMSPPT(IPOROT,IPOINT)
         ENDIF
        ENDDO
       ENDDO
      ENDIF
C
C**** INTERPOLATE FOR MID-SIDE NODES
C
      DO 200 IELEMT=1,NELEMT
       IF(NACTIT.EQ.1) THEN
        IF(LACTIT(IELEMT).EQ.0) GO TO 200     ! skip non active elements
       ENDIF
       LGRUPT=MATNOT(IELEMT)
       LMATST=INT(PROELT(1,LGRUPT))
       LTYPET=INT(PROELT(5,LGRUPT))
       NNODLT=INT(PROELT(2,LGRUPT))
       IF(LTYPET.EQ.101.OR.LTYPET.EQ.104) GO TO 200  !SKIP CONTACT ELEM.
C
C****  CALL THE ELEMENT PROCESSOR TO NOTHING (computes NQUTRT)
C
       CALL ELMLIBT(LNODST(1,IELEMT),PROELT(1,LGRUPT),PROPST(1,LMATST),
     .              ELDATT,ELPRET,ELVART,ELMATT,DUMMYT,WORK1T,    18)
C
       CALL SMONOD(NDIMET,NNODLT,NNODST,NQUTRT)
C
       IF(NNODLT.EQ.NNODST) GO TO 200
C
C**** HEAT FLUXES
C

C
C**** TEMPERATURE GRADIENTS
C

C
C**** INTERNAL VARIABLES
C
       IF(NNUINT.GT.0) THEN
        DO INUINT=1,NNUINT
         DO INODET=1,NNODST
          IPOINT=LNODST(INODET,IELEMT)
          WORK1T(IGSMOT(10)+INODET-1)=SMSTPT(INUINT,IPOINT)
         END DO
         CALL SMOMID(WORK1T(IGSMOT(10)),NDIMET,NNODLT,NQUTRT)
         DO INODET=NNODST+1,NNODLT
          IPOINT=LNODST(INODET,IELEMT)
          SMSTPT(INUINT,IPOINT)=WORK1T(IGSMOT(10)+INODET-1)
         END DO
        END DO
       END IF
C
C**** POROSITY CRITERIA
C
       IF(NPOROT.GT.0.AND.ICALPO.EQ.1) THEN
        DO IPOROT=1,NNUINT
         DO INODET=1,NNODST
          IPOINT=LNODST(INODET,IELEMT)
          WORK1T(IGSMOT(24)+INODET-1)=SMSPPT(IPOROT,IPOINT)
         END DO
         CALL SMOMID(WORK1T(IGSMOT(24)),NDIMET,NNODLT,NQUTRT)
         DO INODET=NNODST+1,NNODLT
          IPOINT=LNODST(INODET,IELEMT)
          SMSPPT(IPOROT,IPOINT)=WORK1T(IGSMOT(24)+INODET-1)
          PREAST(IPOROT,IPOINT)=SMSPPT(IPOROT,IPOINT)
         END DO
        END DO
       END IF
C
  200 CONTINUE
C
      RETURN
      END
