      SUBROUTINE OUTSMOS(ELDATS,ELPRES,ELVARS,ELMATS,LNODSS,MATNOS,
     .                   PROELS,PROPSS,COORDS,PREASS,FPCHAS,DISTOS,
     .                   VELCTS,DISPLS,
     .                   SMSTSS,SMSTNS,ACCPNS,SMSTPS,SMSPPS,
     .                   WORK1S)
C***********************************************************************
C
C**** THIS ROUTINE PERFORMS A SMOOTHING OVER THE VARIABLES GRADIENTS
C
C
C     Notes:
C
C     The call to OUTASS includes NNODES (instead of NNODSS) because the
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
      INCLUDE 'addi_oms.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
      INCLUDE 'nued_om.f'
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
      INCLUDE 'inpo_oms.f'
C
      COMMON/SMOGAUS/FACTAS
      COMMON/POROSI1S/KPOROAS
C
      DIMENSION ACCPNS(NPOINS),        SMSTSS(NSTR1S,NPOINS),
     .          SMSTNS(NSTR1S,NPOINS), SMSTPS(NNUINS,NPOINS),
     .          SMSPPS(NPOROS,NPOINS)
      DIMENSION MATNOS(NELEMS),        LNODSS(NNODES,NELEMS),
     .          PROELS(NPRELS,NGRUPS), PROPSS(NPROPS,NMATSS),
     .          ELDATS(NDATAS),        ELPRES(NPREVS),
     .          ELVARS(NSTATS),        ELMATS(NMATXS)
      DIMENSION COORDS(NDIMES,NPOINS), PREASS(NPOROS,NPOINS),
     .          FPCHAS(NFPCH,NPOINS),  DISTOS(*),
     .          VELCTS(*),             DISPLS(NTOTVMS)
      DIMENSION WORK1S(*)
C
      FACTAS=1.D+00
      KPOROAS=0
C
      DO IPOINS=1,NPOINS
       ACCPNS(IPOINS)=0.D+00
       IF(NMEMO4S.EQ.1) THEN
        CALL RUNENDS('ERROR IN OUTSMOS: NMEMO4S=1 NOT IMPLEMENTED YET')
        DO ISTRES=1,NSTR1S
         SMSTSS(ISTRES,IPOINS)=0.D+00
         SMSTNS(ISTRES,IPOINS)=0.D+00
        ENDDO
       ENDIF
       IF(NNUINS.GT.0) THEN
        CALL RUNENDS('ERROR IN OUTSMOS: NNUINS > 0 NOT IMPLEMENTED YET')
        DO INUINS=1,NNUINS
         SMSTPS(INUINS,IPOINS)=0.D+00
        ENDDO
       ENDIF
       IF(NPOROS.GT.0) THEN
        CALL RUNENDS('ERROR IN OUTSMOS: NPOROS > 0 NOT IMPLEMENTED YET')
        DO IPOROS=1,NPOROS
         SMSPPS(IPOROS,IPOINS)=0.D+00
        ENDDO
       ENDIF
      ENDDO
C
C**** LOOP OVER THE ELEMENTS
C
      DO 100 IELEMS=1,NELEMS
       LGRUPS=MATNOS(IELEMS)
       LMATSS=INT(PROELS(1,LGRUPS))
       LTYPES=INT(PROELS(5,LGRUPS))
       NNODLS=INT(PROELS(2,LGRUPS))
       IF(LTYPES.EQ.101.OR.LTYPES.EQ.104) GO TO 100  !SKIP CONTACT ELEM.
C
C**** READ ELDAT & ELVAR FROM DATA BASE
C
       IF(NMEMO1S.EQ.0.OR.NMEMO2S.EQ.0)
     .  CALL DATBASS(ELDATS,    1,    2)
C
       IF(NMEMO3S.EQ.1.OR.NMEMO4S.EQ.1.OR.NMEMO5S.EQ.0)
     .  CALL DATBASS(ELVARS,    3,    2)    ! current
C
C**** GATHER NODAL COORDINATES
C
       IF(NMEMO1S.EQ.1) THEN    ! coordinates in a global array
        CALL GATHER(COORDS,NDIMES,NPOINS,WORK1S(IGSMOS(19)),NDIMES,
     .              NNODLS,LNODSS(1,IELEMS))
       ENDIF
C
       IF(NMEMO4S.EQ.1) THEN          ! not used now !
        IF(ITERME.GT.0) THEN          ! bidirectional coupled
         IF(ITERMD.EQ.1) THEN         ! deformed shape
C
C**** SCALAR GATHER OPERATIONS ( DISPLT ---> WORK1T(IGSMOT(30)) )
C
          CALL GATHER(DISPLS,NDOFCMS,NPOINS,WORK1S(IGSMOS(30)),NDOFCMS,
     .                NNODLS,
     .                LNODSS(1,IELEMS))
         ENDIF                        ! itermd.eq.1
        ENDIF                         ! iterme.gt.0
C
C**** CALL THE ELEMENT PROCESSOR TO EVALUATE ELEMENT CONTRIBUTION TO
C     GLOBAL VECTOR OF HEAT FLUXES
C
        CALL ELMLIBS(LNODSS(1,IELEMS),PROELS(1,LGRUPS),PROPSS(1,LMATSS),
     .               ELDATS,ELPRES,ELVARS,ELMATS,WORK1S,      12)
C
C**** COMPUTES ACCPNS
C
        DO INODES=1,NNODSS
         LPOINS=LNODSS(INODES,IELEMS)
         IF(LPOINS.NE.0) ACCPNS(LPOINS)=ACCPNS(LPOINS)+FACTAS
        ENDDO
C
C**** INCLUDE IT IN A GLOBAL VECTOR
C
        CALL OUTASS(LNODSS,SMSTSS,WORK1S(IGSMOS(5)),FACTAS,
     .              NNODES,IELEMS,NELEMS,NPOINS,
     .              NSTR1S)
C
C**** CALL THE ELEMENT PROCESSOR TO EVALUATE ELEMENT CONTRIBUTION TO
C     GLOBAL VECTOR OF VARIABLES GRADIENTS
C
        CALL ELMLIBS(LNODSS(1,IELEMS),PROELS(1,LGRUPS),PROPSS(1,LMATSS),
     .               ELDATS,ELPRES,ELVARS,ELMATS,WORK1S,      13)
C
C**** INCLUDE IT IN A GLOBAL VECTOR
C
        CALL OUTASS(LNODSS,SMSTNS,WORK1S(IGSMOS(7)),FACTAS,
     .              NNODES,IELEMS,NELEMS,NPOINS,
     .              NSTR1S)
       ENDIF     ! nmemo4.eq.1
C
c      IF(NNUINT.GT.0) THEN
c       IF(NMEMO4.EQ.0) THEN
c        IF(ITERME.GT.0) THEN          ! bidirectional coupled
c         IF(ITERMD.EQ.1) THEN         ! deformed shape
C
C**** SCALAR GATHER OPERATIONS ( DISPLT ---> WORK1T(IGSMOT(30)) )
C
c          CALL GATHER(DISPLT,NDOFCM,NPOINT,WORK1T(IGSMOT(30)),NDOFCM,
c    .                 NNODLT,
c    .                 LNODST(1,IELEMT))
c         ENDIF                        ! itermd.eq.1
c        ENDIF                         ! iterme.gt.0
c       ENDIF                          ! nmemo4.eq.0
C
C**** CALL THE ELEMENT PROCESSOR TO EVALUATE ELEMENT CONTRIBUTION TO
C     GLOBAL VECTOR OF INTERNAL VARIABLES
C
c       CALL ELMLIBT(LNODST(1,IELEMT),PROELT(1,LGRUPT),PROPST(1,LMATST),
c    .               ELDATT,ELPRET,ELVART,ELMATT,WORK1T,   14)
C
C**** COMPUTES ACCPN (IF NECESSARY)
C
c       IF(NMEMO4.EQ.0) THEN
c        DO INODET=1,NNODST
c         LPOINT=LNODST(INODET,IELEMT)
c         IF(LPOINT.NE.0) ACCPNT(LPOINT)=ACCPNT(LPOINT)+FACTAT
c        ENDDO
c       ENDIF
C
C**** INCLUDE IT IN A GLOBAL VECTOR
C
c       CALL OUTASS(LNODST,SMSTPT,WORK1T(IGSMOT(10)),FACTAT,
c    .              NNODET,IELEMT,NELEMT,NPOINT,
c    .              NNUINT)
C
c      ENDIF     ! nnuint.gt.0
C
C**** EVALUATE ELEM. CONTRIBUTION TO GLOBAL VECTOR OF POROSITY CRITERIA
C
c      IF(NPOROT.GT.0) THEN
c       IF(IMICR.EQ.0) THEN
C
C**** GATHER NODAL POROSITY CRITERIA
C
c        NGAULT=INT(PROELT(4,LGRUPT))
c        NNODGT=NNODLT
c        IF(NGAULT.GT.NNODLT) NNODGT=NGAULT
c        CALL GATHER(PREAST,NPOROT,NPOINT,WORK1T(IGSMOT(25)),NPOROT,
c    .               NNODGT,LNODST(1,IELEMT))
C
c        IF(NMEMO5.EQ.0) THEN
C
C**** SCALAR GATHER OPERATIONS ( DISTOT ---> ELDIST(ELVART) )
C
c         CALL GATHER(DISTOT,NDOFCT,NPOINT,ELVART,NDOFCT,NNODLT,
c    .                LNODST(1,IELEMT))
c        ELSE
C
C**** SCALAR GATHER OPERATIONS ( DISTOT ---> ELDIST )
C
c         CALL GATHER(DISTOT,NDOFCT,NPOINT,WORK1T(IGSMOT(21)),NDOFCT,
c    .                NNODLT,
c    .                LNODST(1,IELEMT))
c        ENDIF                   ! nmemo5.eq.0
C
C**** GATHER NODAL PHASE-CHANGE FUNCTION (MACROSCOPIC)
C
c        CALL GATHER(FPCHAT,NFPCH,NPOINT,WORK1T(IGSMOT(26)),NFPCH,
c    .               NNODLT,LNODST(1,IELEMT))
C
C**** SCALAR GATHER OPERATIONS ( VELCTT ---> VELCMT ) - TEMP. RATE
C
c        IF(KDYNAT.EQ.1)
c    .    CALL GATHER(VELCTT,NDOFCT,NPOINT,WORK1T(IGSMOT(27)),NDOFCT,
c    .                NNODLT,
c    .                LNODST(1,IELEMT))
C
c        IF(NMEMO4.EQ.0.AND.NNUINT.EQ.0) THEN
c         IF(ITERME.GT.0) THEN          ! bidirectional coupled
c          IF(ITERMD.EQ.1) THEN         ! deformed shape
C
C**** SCALAR GATHER OPERATIONS ( DISPLT ---> WORK1T(IGSMOT(30)) )
C
c           CALL GATHER(DISPLT,NDOFCM,NPOINT,WORK1T(IGSMOT(30)),NDOFCM,
c    .                  NNODLT,
c    .                  LNODST(1,IELEMT))
c          ENDIF                        ! itermd.eq.1
c         ENDIF                         ! iterme.gt.0
c        ENDIF                          ! nmemo4.eq.0.and.nnuint.eq.0
C
C**** CALL THE ELEMENT PROCESSOR TO EVALUATE ELEMENT CONTRIBUTION TO
C     GLOBAL VECTOR OF POROSITY CRITERIA
C
c        CALL ELMLIBT(LNODST(1,IELEMT),PROELT(1,LGRUPT),
c    .                PROPST(1,LMATST),
c    .                ELDATT,ELPRET,ELVART,ELMATT,WORK1T,   16)
C
C**** ADDITIONAL CONTROL
C
c        KPOROA=KPOROA+KPOROT
C
C**** COMPUTES ACCPN (IF NECESSARY)
C
c        IF(NMEMO4.EQ.0.AND.NNUINT.EQ.0) THEN
c         DO INODET=1,NNODST
c          LPOINT=LNODST(INODET,IELEMT)
c          IF(LPOINT.NE.0) ACCPNT(LPOINT)=ACCPNT(LPOINT)+FACTAT
c         ENDDO
c        ENDIF
C
C**** INCLUDE IT IN A GLOBAL VECTOR
C
c        CALL OUTASS(LNODST,SMSPPT,WORK1T(IGSMOT(24)),FACTAT,
c    .               NNODET,IELEMT,NELEMT,NPOINT,
c    .               NPOROT)
C
c       ENDIF    ! imicr.eq.0
c      ENDIF     ! nporot.gt.0
C
  100 CONTINUE
C
C**** ASSIGNS MICROSCOPICAL PHASE-CHANGE FUNCTION TO FPCHAT
C
C     Notes:
C
C     IPLUAT(INUPM) (defined in pointes.f and passed through auxl_oms.f)
C     gives the position of the phase-change functions in
C     SMSTPT(IPLUAT(INUPM),IPOINT)
C
C     FPCHAT can be transferred to mechanical problem (coupled problems)
C     or be used in the porosity criteria (uncouped & coupled problems)
C
c     IF(IMICR.EQ.1) THEN
c      IF(NNUPM.GT.0) THEN
c       DO INUPM=1,NNUPM
c        INUAU=IPLUAT(INUPM)
c        DO IPOINT=1,NPOINT
c         FAUXX=FPCHAT(INUPM+NNUPC,IPOINT)              ! last converged
c         IF(ACCPNT(IPOINT).NE.0.0)
c    .     FPCHAT(INUPM+NNUPC,IPOINT)=SMSTPT(INUAU,IPOINT)/    ! current
c    .                                                    ACCPNT(IPOINT)
c         IF(IPRCOT.GT.0)
c    .     FPCHAT(INUPM+NNUPC+NNUPT,IPOINT)=
c    .     (FPCHAT(INUPM+NNUPC,IPOINT)-FAUXX)/DTIMET          ! dot f_pc
c        ENDDO
c       ENDDO
c      ENDIF            ! nnupm.gt.0
C
C**** ASSIGNS OTHER VARIABLES TO FPCHAT
C
C     Notes:
C
C     IPLUOT(INUPO) (defined in pointes.f and passed through auxl_oms.f)
C     gives the position of the other variables in
C     SMSTPT(IPLUOT(INUPO),IPOINT)
C
C     FPCHAT can be transferred to mechanical problem (coupled problems)
C     or be used in the porosity criteria (uncouped & coupled problems)
C
c      IF(NNUPO.GT.0) THEN
c       DO INUPO=1,NNUPO
c        INUAU=IPLUOT(INUPO)
c        DO IPOINT=1,NPOINT
c         IF(ACCPNT(IPOINT).NE.0.0)
c    .     FPCHAT(INUPO+2*NNUPT,IPOINT)=SMSTPT(INUAU,IPOINT)/
c    .                                                    ACCPNT(IPOINT)
c        ENDDO
c       ENDDO
c      ENDIF            ! nnupo.gt.0
C
C**** LOOP OVER THE ELEMENTS
C
c      DO 101 IELEMT=1,NELEMT
c       LGRUPT=MATNOT(IELEMT)
c       LMATST=INT(PROELT(1,LGRUPT))
c       LTYPET=INT(PROELT(5,LGRUPT))
c       NNODLT=INT(PROELT(2,LGRUPT))
c       IF(LTYPET.EQ.101.OR.LTYPET.EQ.104) GO TO 101 !SKIP CONTACT ELEM.
C
C**** READ ELDAT & ELVAR FROM DATA BASE
C
c       IF(NMEMO1.EQ.0.OR.NMEMO2.EQ.0)
c    .   CALL DATBAST(ELDATT,    1,    2)
C
c       IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1)
c    .   CALL DATBAST(ELVART,    3,    2)    ! current
C
C**** GATHER NODAL COORDINATES
C
c       IF(NMEMO1.EQ.1) THEN    ! coordinates in a global array
c        CALL GATHER(COORDT,NDIMET,NPOINT,WORK1T(IGSMOT(19)),NDIMET,
c    .               NNODLT,LNODST(1,IELEMT))
c       ENDIF
C
C**** EVALUATE ELEM. CONTRIBUTION TO GLOBAL VECTOR OF POROSITY CRITERIA
C
C     Note: the porosity criteria may contain microstructural
C           phase-change functions and other microstructural variables
C
c       IF(NPOROT.GT.0) THEN
C
C**** GATHER NODAL POROSITY CRITERIA
C
c        NGAULT=INT(PROELT(4,LGRUPT))
c        NNODGT=NNODLT
c        IF(NGAULT.GT.NNODLT) NNODGT=NGAULT
c        CALL GATHER(PREAST,NPOROT,NPOINT,WORK1T(IGSMOT(25)),NPOROT,
c    .               NNODGT,LNODST(1,IELEMT))
C
C**** GATHER NODAL PHASE-CHANGE FUNCTION (MACROSCOPIC, MICROSCOPIC &
C     OTHER MICROSTRUCTURE VARIABLES)
C
c        CALL GATHER(FPCHAT,NFPCH,NPOINT,WORK1T(IGSMOT(26)),NFPCH,
c    .               NNODLT,LNODST(1,IELEMT))
C
c        IF(NMEMO4.EQ.0.AND.NNUINT.EQ.0) THEN
c         IF(ITERME.GT.0) THEN          ! bidirectional coupled
c          IF(ITERMD.EQ.1) THEN         ! deformed shape
C
C**** SCALAR GATHER OPERATIONS ( DISPLT ---> WORK1T(IGSMOT(30)) )
C
c           CALL GATHER(DISPLT,NDOFCM,NPOINT,WORK1T(IGSMOT(30)),NDOFCM,
c    .                  NNODLT,
c    .                  LNODST(1,IELEMT))
c          ENDIF                        ! itermd.eq.1
c         ENDIF                         ! iterme.gt.0
c        ENDIF                          ! nmemo4.eq.0.and.nnuint.eq.0
C
C**** CALL THE ELEMENT PROCESSOR TO EVALUATE ELEMENT CONTRIBUTION TO
C     GLOBAL VECTOR OF POROSITY CRITERIA
C
c        CALL ELMLIBT(LNODST(1,IELEMT),PROELT(1,LGRUPT),
c    .                PROPST(1,LMATST),
c    .                ELDATT,ELPRET,ELVART,ELMATT,WORK1T,   16)
C
C**** ADDITIONAL CONTROL
C
c        KPOROA=KPOROA+KPOROT
C
C**** COMPUTES ACCPN (IF NECESSARY)
C
c        IF(NMEMO4.EQ.0.AND.NNUINT.EQ.0) THEN
c         DO INODET=1,NNODST
c          LPOINT=LNODST(INODET,IELEMT)
c          IF(LPOINT.NE.0) ACCPNT(LPOINT)=ACCPNT(LPOINT)+FACTAT
c         ENDDO
c        ENDIF
C
C**** INCLUDE IT IN A GLOBAL VECTOR
C
c        CALL OUTASS(LNODST,SMSPPT,WORK1T(IGSMOT(24)),FACTAT,
c    .               NNODET,IELEMT,NELEMT,NPOINT,
c    .               NPOROT)
C
c       ENDIF     ! nporot.gt.0
C
c 101  CONTINUE
c     ENDIF             ! imicr.eq.1
C
C**** PERFORMS AVERAGE MEAN FOR HEAT FLUXES
C
c     IF(NMEMO4.EQ.1) THEN
c      DO IPOINT=1,NPOINT
c       DO ISTRET=1,NSTR1T
c        IF(ACCPNT(IPOINT).NE.0.0)
c    .    SMSTST(ISTRET,IPOINT)=SMSTST(ISTRET,IPOINT)/ACCPNT(IPOINT)
c       ENDDO
c      ENDDO
c     ENDIF
C
C**** PERFORMS AVERAGE MEAN FOR TEMPERATURE GRADIENTS (not implem. yet)
C

C
C**** PERFORMS AVERAGE MEAN FOR INTERNAL VARIABLES
C
c     IF(NNUINT.GT.0) THEN
c      DO IPOINT=1,NPOINT
c       DO INUINT=1,NNUINT
c        IF(ACCPNT(IPOINT).NE.0.0)
c    .    SMSTPT(INUINT,IPOINT)=SMSTPT(INUINT,IPOINT)/ACCPNT(IPOINT)
c       ENDDO
c      ENDDO
c     ENDIF
C
C**** PERFORMS AVERAGE MEAN FOR POROSITY CRITERIA & ASSIGNS TO PREAST
C
c     IF(NPOROT.GT.0) THEN
c      DO IPOINT=1,NPOINT
c       DO IPOROT=1,NPOROT
c        IF(ACCPNT(IPOINT).NE.0.0) THEN
c         SMSPPT(IPOROT,IPOINT)=SMSPPT(IPOROT,IPOINT)/ACCPNT(IPOINT)
c         PREAST(IPOROT,IPOINT)=SMSPPT(IPOROT,IPOINT)
c        ENDIF
c       ENDDO
c      ENDDO
c     ENDIF
C
C**** INTERPOLATE FOR MID-SIDE NODES
C
c     DO 200 IELEMT=1,NELEMT
c      LGRUPT=MATNOT(IELEMT)
c      LMATST=INT(PROELT(1,LGRUPT))
c      LTYPET=INT(PROELT(5,LGRUPT))
c      NNODLT=INT(PROELT(2,LGRUPT))
c      IF(LTYPET.EQ.101.OR.LTYPET.EQ.104) GO TO 200  !SKIP CONTACT ELEM.
C
C**** CALL THE ELEMENT PROCESSOR TO "NOTHING" (computes NQUTRS)
C
c      CALL ELMLIBT(LNODST(1,IELEMT),PROELT(1,LGRUPT),
c    .              PROPST(1,LMATST),
c    .              ELDATT,ELPRET,ELVART,ELMATT,WORK1T,   17)
C
c      CALL SMONOD(NDIMET,NNODLT,NNODST,NQUTRS)
C
c      IF(NNODLT.EQ.NNODST) GO TO 200
C
C**** HEAT FLUXES
C

C
C**** TEMPERATURE GRADIENTS
C

C
C**** INTERNAL VARIABLES
C
c      IF(NNUINT.GT.0) THEN
c       DO INUINT=1,NNUINT
c        DO INODET=1,NNODST
c         IPOINT=LNODST(INODET,IELEMT)
c         WORK1T(IGSMOT(10)+INODET-1)=SMSTPT(INUINT,IPOINT)
c        END DO
c        CALL SMOMID(WORK1T(IGSMOT(10)),NDIMET,NNODLT,NQUTRS)
c        DO INODET=NNODST+1,NNODLT
c         IPOINT=LNODST(INODET,IELEMT)
c         SMSTPT(INUINT,IPOINT)=WORK1T(IGSMOT(10)+INODET-1)
c        END DO
c       END DO
c      END IF
C
C**** POROSITY CRITERIA
C
c      IF(NPOROT.GT.0) THEN
c       DO IPOROT=1,NNUINT
c        DO INODET=1,NNODST
c         IPOINT=LNODST(INODET,IELEMT)
c         WORK1T(IGSMOT(24)+INODET-1)=SMSPPT(IPOROT,IPOINT)
c        END DO
c        CALL SMOMID(WORK1T(IGSMOT(24)),NDIMET,NNODLT,NQUTRS)
c        DO INODET=NNODST+1,NNODLT
c         IPOINT=LNODST(INODET,IELEMT)
c         SMSPPT(IPOROT,IPOINT)=WORK1T(IGSMOT(24)+INODET-1)
c        END DO
c       END DO
c      END IF
C
c 200 CONTINUE
C
      RETURN
      END
