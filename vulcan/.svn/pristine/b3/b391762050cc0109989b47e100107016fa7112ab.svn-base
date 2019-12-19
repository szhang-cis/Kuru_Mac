      SUBROUTINE SETMTXT(ELDATT,ELPRET,ELVART,ELMATT,LNODST,MATNOT,
     .                   PROELT,PROPST,COORDT,TEMPIT,FPCHAT,DISPLT,
     .                   WORK1T)
C***********************************************************************
C
C**** THIS ROUTINE SETS UP SOME CONSTANT MATRICES FOR FUTURE USE:
C
C         - ROTATION MATRIX             (RMAT1)
C         - ELASTIC CONSTANTS MATRIX    (EPMTX)
C         - NODAL SMOOTHING MATRICES    (SMATR;RMATR,FMATR)
C         - GAUSSIAN SMOOTHING MATRICES (EMASS)
C
C     IT ALSO EVALUATES ONCE AND FOR ALL:
C
C         - SHAPE FUNCTIONS          (SHAPE)
C         - CARTESIAN DERIVATIVES    (CARTD)
C         - DIFFERENTIAL VOLUMEN     (DVOLU)
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
      COMMON/JACOBSTA/IERORT,KERORT
C
      DIMENSION MATNOT(NELEMT),        LNODST(NNODET,NELEMT),
     .          PROELT(NPRELT,NGRUPT), PROPST(NPROPT,NMATST),
     .          ELDATT(NDATAT),        ELPRET(NPREVT),
     .          ELVART(NSTATT),        ELMATT(NMATXT),
     .          COORDT(NDIMET,NPOINT), TEMPIT(NPOINT),
     .          FPCHAT(NFPCH,NPOINT),  DISPLT(NTOTVM),
     .          WORK1T(*)
C
      CALL CPUTIMT(TIME1T)
C
C**** LOOP ON ELEMENTS
C
C     Note: KELEMT is used instead of IELEMT as a loop index in order to
C           avoid compiling problems when NOCOIT>0 is considered since
C           in this case a JELEMT do loop is also performed and IELEMT
C           is reserved to be used in datbast.f for these two loops
C           (IELEMT is transferred by auxl_omt.f).
C
      KERORT=0
      DO 1000 KELEMT=1,NELEMT
      IELEMT=KELEMT
      LGRUPT=MATNOT(KELEMT)
      LMATST=INT(PROELT(1,LGRUPT))
      NNODLT=INT(PROELT(2,LGRUPT))
      LTYPET=INT(PROELT(5,LGRUPT))
C
      IERORT=0
C
      IF(NMEMO1.EQ.0) THEN     ! coordinates in an elemental array
C
C**** READ ELDAT ( coordinates ) FROM DATA-BASE
C
       CALL DATBAST(ELDATT,    1,    2)
C
      ELSE                     ! coordinates in a global array
       CALL GATHER(COORDT,NDIMET,NPOINT,WORK1T(ISETMT(6)),NDIMET,
     .             NNODLT,LNODST(1,KELEMT))
      ENDIF
C
C**** READ ELVAR FOR MICROSTRUCTURAL PROBLEMS
C
      IF(IMICR.EQ.1) THEN
       CALL DATBAST(ELVART,    3,    2)
      ENDIF
C
C**** SCALAR GATHER OPERATIONS ( TEMPIT ---> WORK1T(ISETMT(13)) )
C
      CALL GATHER(TEMPIT,NDOFCT,NPOINT,WORK1T(ISETMT(13)),NDOFCT,NNODLT,
     .            LNODST(1,KELEMT))     ! initial temperature
C
      IF(ITERME.GT.0) THEN          ! bidirectional coupled
       IF(ITERMD.EQ.1) THEN         ! deformed shape
C
C**** SCALAR GATHER OPERATIONS ( DISPLT ---> WORK1T(ISETMT(18)) )
C
        CALL GATHER(DISPLT,NDOFCM,NPOINT,WORK1T(ISETMT(18)),NDOFCM,
     .              NNODLT,
     .              LNODST(1,KELEMT))
       ENDIF                        ! itermd.eq.1
      ENDIF                         ! iterme.gt.0
C
C**** GATHER PHASE-CHANGE VECTOR ( FPCHAT ---> WORK1T(ISETMT(14)) )
C
      IF(LTYPET.EQ.5)
     . CALL GATHER(FPCHAT,NFPCH,NPOINT,WORK1T(ISETMT(14)),NFPCH,
     .             NNODLT,LNODST(1,KELEMT))
C
C**** SET UP ARRAY ELDAT & WRITE IT TO DATA BASE
C
      CALL ELMLIBT(LNODST(1,KELEMT),PROELT(1,LGRUPT),PROPST(1,LMATST),
     .             ELDATT,ELPRET,ELVART,ELMATT,DUMMYT,WORK1T,     1)
C
      IF(IERORT.NE.0) GO TO 1000
C
C**** WRITE ELDAT TO DATA-BASE
C
      IF(NMEMO1.EQ.0.OR.NMEMO2.EQ.0)
     . CALL DATBAST(ELDATT,    1,    1)
C
C**** WRITE ELVAR FOR MICROSTRUCTURAL PROBLEMS (initial conditions)
C
C     ELVAR will be considered in adddatt.f as last converged values
C
C     ELVAR IS NOT WRITTEN FOR RESTART PROBLEMS (the initial conditions
C     are in *.pan file)
C
      IF(INITIT.EQ.1) THEN
       IF(IMICR.EQ.1)
     .  CALL DATBAST(ELVART,    3,    1)
      ENDIF
C
      IF(LTYPET.EQ.5)
     . CALL SCATERA(WORK1T(ISETMT(14)),NFPCH,NNODLT,FPCHAT,NFPCH,
     .              NPOINT,LNODST(1,KELEMT))
C
C**** DEALS WITH CONTACT ELEMENTS FOR NON-COINCIDENT MESH
C
      IF(NOCOIT.GT.0) THEN
       IF(LTYPET.NE.104) GO TO 1000          ! skip non-contact elements
       NOCOLT=INT(PROELT(11,LGRUPT))
       IF(NOCOLT.EQ.0) GO TO 1000            ! skip coincident cont. el.
C
       NGAULT=INT(PROELT(4,LGRUPT))
C
       ISETCT=INT(PROELT(8,LGRUPT))          ! contact set
       IMAEST=INT(PROELT(9,LGRUPT))          ! contact index
C
C**** RECOVERS COORDINATES & NORMAL OF MASTER CONTACT ELEMENT (KELEMT)
C     AT GAUSS POINTS
C
       IF(NMEMO1.EQ.0) THEN
        IAUXXT=-1
        DO IGAUST=1,NGAULT
         DO IDIMET=1,NDIMET
          IAUXXT=IAUXXT+1
          WORK1T(ISETMT(19)+IAUXXT)=ELDATT(IDATAT( 4)+IAUXXT)    ! GPCOD
          WORK1T(ISETMT(21)+IAUXXT)=ELDATT(IDATAT(10)+IAUXXT)    ! VNORL
         ENDDO
        ENDDO
       ELSE
        IAUXXT=-1
        DO IGAUST=1,NGAULT
         DO IDIMET=1,NDIMET
          IAUXXT=IAUXXT+1
          WORK1T(ISETMT(19)+IAUXXT)=WORK1T(ISETMT( 9)+IAUXXT)    ! GPCOD
          WORK1T(ISETMT(21)+IAUXXT)=WORK1T(ISETMT(15)+IAUXXT)    ! VNORL
         ENDDO
        ENDDO
       ENDIF
C
       DO IGAUST=1,NGAULT             ! proper initialization for solver
        DO INODLT=1,NNODLT
         INODXT=NNODLT*IGAUST+INODLT
         LNODST(NNODLT*IGAUST+INODLT,KELEMT)=LNODST(INODLT,KELEMT)
        ENDDO
       ENDDO
C
       DO IGAUST=1,NGAULT
        WORK1T(ISETMT(22)+IGAUST-1)=-1.0D+10 ! GAPMI: min. dist. along n
       ENDDO
C
C**** LOOP ON ELEMENTS
C
C     Note: this loop is only necessary to compute nodal coordinates of
C           the slave elements
C
       DO 2000 JELEMT=1,NELEMT
       IF(KELEMT.EQ.JELEMT) GO TO 2000       ! skip equal elements
       IELEMT=JELEMT
       JGRUPT=MATNOT(JELEMT)
       JMATST=INT(PROELT(1,JGRUPT))
       JNODLT=INT(PROELT(2,JGRUPT))
       JTYPET=INT(PROELT(5,JGRUPT))
       IF(JTYPET.NE.104) GO TO 2000          ! skip non-contact elements
       JOCOLT=INT(PROELT(11,JGRUPT))
       IF(JOCOLT.EQ.0) GO TO 2000            ! skip coincident cont. el.
C
       JSETCT=INT(PROELT(8,JGRUPT))          ! contact set
       JMAEST=INT(PROELT(9,JGRUPT))          ! contact index
C
C**** CONTACT CRITERIA
C
       IAUXXT=0
       IF((IMAEST.EQ.1.AND.JMAEST.EQ.1).OR.  ! no self-contact (nsc)
     .    (IMAEST.EQ.1.AND.JMAEST.EQ.3).OR.
     .    (IMAEST.EQ.3.AND.JMAEST.EQ.1)) THEN
        IF(ISETCT.EQ.JSETCT) THEN
         IF(LGRUPT.NE.JGRUPT) IAUXXT=1
        ENDIF
       ENDIF
       IF((IMAEST.EQ.2.AND.JMAEST.EQ.2).OR.  ! self-contact (sc)
     .    (IMAEST.EQ.2.AND.JMAEST.EQ.3).OR.
     .    (IMAEST.EQ.3.AND.JMAEST.EQ.2)) THEN
        IF(ISETCT.EQ.JSETCT) THEN
         IF(LGRUPT.EQ.JGRUPT) IAUXXT=1
        ENDIF
       ENDIF
       IF(IMAEST.EQ.3.AND.JMAEST.EQ.3) THEN  ! both (nsc & sc)
        IF(ISETCT.EQ.JSETCT) IAUXXT=1
       ENDIF
       IF(IAUXXT.EQ.0) GO TO 2000         ! skip due to contact criteria
C
       IF(NMEMO1.EQ.0) THEN     ! coordinates in an elemental array
C
C**** READ ELDAT ( coordinates ) FROM DATA-BASE
C
        CALL DATBAST(ELDATT,    1,    2)
C
       ELSE                     ! coordinates in a global array
        CALL GATHER(COORDT,NDIMET,NPOINT,WORK1T(ISETMT(6)),NDIMET,
     .              NNODLT,LNODST(1,JELEMT))
       ENDIF
C
       IF(ITERME.GT.0) THEN          ! bidirectional coupled
        IF(ITERMD.EQ.1) THEN         ! deformed shape
C
C**** SCALAR GATHER OPERATIONS ( DISPLT ---> WORK1T(ISETMT(18)) )
C
         CALL GATHER(DISPLT,NDOFCM,NPOINT,WORK1T(ISETMT(18)),NDOFCM,
     .               NNODLT,
     .               LNODST(1,JELEMT))
        ENDIF                        ! itermd.eq.1
       ENDIF                         ! iterme.gt.0
C
C**** SET UP ARRAY ELDAT & WRITE IT TO DATA BASE
C
       CALL ELMLIBT(LNODST(1,JELEMT),PROELT(1,JGRUPT),PROPST(1,JMATST),
     .              ELDATT,ELPRET,ELVART,ELMATT,DUMMYT,WORK1T,     1)
C
C**** RECOVERS COORDINATES OF SLAVE CONTACT ELEMENT (JELEMT) AT NODES
C
       IF(NMEMO1.EQ.0) THEN
        IAUXXT=-1
        DO IGAUST=1,NGAULT
         DO IDIMET=1,NDIMET
          IAUXXT=IAUXXT+1
          WORK1T(ISETMT(20)+IAUXXT)=ELDATT(IDATAT( 1)+IAUXXT)    ! ELCOD
         ENDDO
        ENDDO
       ELSE
        IAUXXT=-1
        DO INODLT=1,NNODLT
         DO IDIMET=1,NDIMET
          IAUXXT=IAUXXT+1
          WORK1T(ISETMT(20)+IAUXXT)=WORK1T(ISETMT( 6)+IAUXXT)    ! ELCOD
         ENDDO
        ENDDO
       ENDIF
C
C**** TASKS (ACCORDING TO NDIMET):
C
C     1) COMPUTES COORDINATES OF SLAVE CONTACT POINT CORRESPONDING TO
C        THE GAUSS POINT OF THE MASTER ELEMENT (THESE TWO POINTS ARE
C        THOSE PRESENTING THE MINIMUN DISTANCE BETWEEN THEM)
C
C     2) REDEFINES CONNECTIVITY ARRAY OF THE MASTER ELEMENT BY ADDING
C        NODES OF THE SLAVE ELEMEMT
C
C     3) COMPUTES NATURAL COORDINATES OF CONTACT POINT & ASSIGNS THEM TO
C        ELVART (STRSGT)
C
C
C        Notes:
C
C        Index IAUXZT of STRSGT(NDIMET-1): natural coordinates of
C        contact point
C
C        The condition NMEMO4.EQ.1 needed to call datbast.f with ELVART
C        has been checked in proinpt.f
C
       IF(NDIMET.EQ.2) THEN
        X1A=WORK1T(ISETMT(20)  )
        Y1A=WORK1T(ISETMT(20)+1)
        X1B=WORK1T(ISETMT(20)+2)
        Y1B=WORK1T(ISETMT(20)+3)
        DO IGAUST=1,NGAULT
         X1= WORK1T(ISETMT(19)+(IGAUST-1)*NDIMET  )
         Y1= WORK1T(ISETMT(19)+(IGAUST-1)*NDIMET+1)
         PX1=WORK1T(ISETMT(21)+(IGAUST-1)*NDIMET  )
         PY1=WORK1T(ISETMT(21)+(IGAUST-1)*NDIMET+1)
C
         CALL CPCALC(X1A,Y1A,Z1A,X1B,Y1B,Z1B,X1C,Y1C,Z1C,X1D,Y1D,Z1D,
     .               X1, Y1, Z1, PX1,PY1,PZ1,
     .               NDIMET,NNODLT,RCOOR,SCOOR,DISTT,DISMI,IAUXET,
     .               RCINF,SCINF,RCSUP)
         IF(IAUXET.EQ.1)
     .    CALL RUNENDT('ERROR: NATURAL COORD. NOT FOUND IN cpcalc.f')
C
         GAPMI=WORK1T(ISETMT(22)+IGAUST-1)
         IF(RCOOR.GE.-1.0D0.AND.RCOOR.LE.1.0D0.AND.
     .                                            (DISMI.GT.GAPMI)) THEN
          WORK1T(ISETMT(22)+IGAUST-1)=DISMI
          IAUXYT=(IGAUST-1)*NNODLT
          LNODST(NNODLT+IAUXYT+1,KELEMT)=LNODST(1,JELEMT)
          LNODST(NNODLT+IAUXYT+2,KELEMT)=LNODST(2,JELEMT)
C
          IELEMT=KELEMT
          IF(IPRCOT.EQ.0) THEN
           CALL DATBAST(ELVART,    3,    2)    ! READ CURRENT VALUES
          ELSE
           CALL DATBAST(ELVART,    5,    2)    ! READ LAST CONV. VALUES
          ENDIF
C
          IAUXZT=(IGAUST-1)*NSTR1T
          ELVART(ISTATT(4)+IAUXZT)=RCOOR
C
          IF(IPRCOT.EQ.0) THEN
           CALL DATBAST(ELVART,    3,    1)    ! WRITE CURRENT VALUES
          ELSE
           CALL DATBAST(ELVART,    5,    1)    ! WRITE LAST CONV. VALUES
          ENDIF
          IELEMT=JELEMT
         ENDIF            ! rcoor.ge.-1.0d0....
C
        ENDDO             ! igaust=1,ngault
       ENDIF              ! ndimet.eq.2
C
       IF(NDIMET.EQ.3) THEN
        X1A=WORK1T(ISETMT(20)  )
        Y1A=WORK1T(ISETMT(20)+1)
        Z1A=WORK1T(ISETMT(20)+2)
        X1B=WORK1T(ISETMT(20)+3)
        Y1B=WORK1T(ISETMT(20)+4)
        Z1B=WORK1T(ISETMT(20)+5)
        X1C=WORK1T(ISETMT(20)+6)
        Y1C=WORK1T(ISETMT(20)+7)
        Z1C=WORK1T(ISETMT(20)+8)
        IF(NNODLT.GT.3) THEN
         X1D=WORK1T(ISETMT(20)+ 9)
         Y1D=WORK1T(ISETMT(20)+10)
         Z1D=WORK1T(ISETMT(20)+11)
        ENDIF
        DO IGAUST=1,NGAULT
         X1= WORK1T(ISETMT(19)+(IGAUST-1)*NDIMET  )
         Y1= WORK1T(ISETMT(19)+(IGAUST-1)*NDIMET+1)
         Z1= WORK1T(ISETMT(19)+(IGAUST-1)*NDIMET+2)
         PX1=WORK1T(ISETMT(21)+(IGAUST-1)*NDIMET  )
         PY1=WORK1T(ISETMT(21)+(IGAUST-1)*NDIMET+1)
         PZ1=WORK1T(ISETMT(21)+(IGAUST-1)*NDIMET+2)
C
         CALL CPCALC(X1A,Y1A,Z1A,X1B,Y1B,Z1B,X1C,Y1C,Z1C,X1D,Y1D,Z1D,
     .               X1, Y1, Z1, PX1,PY1,PZ1,
     .               NDIMET,NNODLT,RCOOR,SCOOR,DISTT,DISMI,IAUXET,
     .               RCINF,SCINF,RCSUP)
         IF(IAUXET.EQ.1)
     .    CALL RUNENDT('ERROR: NATURAL COORD. NOT FOUND IN cpcalc.f')
C
         GAPMI=WORK1T(ISETMT(22)+IGAUST-1)
         IF((RCOOR.GE.RCINF.AND.RCOOR.LE.1.0D0).AND.
     .      (SCOOR.GE.SCINF.AND.SCOOR.LE.(1.0D0-RCSUP)).AND.
     .                                            (DISMI.GT.GAPMI)) THEN
          WORK1T(ISETMT(22)+IGAUST-1)=DISMI
          IAUXYT=(IGAUST-1)*NNODLT
          LNODST(NNODLT+IAUXYT+1,KELEMT)=LNODST(1,JELEMT)
          LNODST(NNODLT+IAUXYT+2,KELEMT)=LNODST(2,JELEMT)
          LNODST(NNODLT+IAUXYT+3,KELEMT)=LNODST(3,JELEMT)
          IF(NNODLT.GT.3)
     .    LNODST(NNODLT+IAUXYT+4,KELEMT)=LNODST(4,JELEMT)
C
          IELEMT=KELEMT
          IF(IPRCOT.EQ.0) THEN
           CALL DATBAST(ELVART,    3,    2)    ! READ CURRENT VALUES
          ELSE
           CALL DATBAST(ELVART,    5,    2)    ! READ LAST CONV. VALUES
          ENDIF
C
          IAUXZT=(IGAUST-1)*NSTR1T
          ELVART(ISTATT(4)+IAUXZT  )=RCOOR
          ELVART(ISTATT(4)+IAUXZT+1)=SCOOR
C
          IF(IPRCOT.EQ.0) THEN
           CALL DATBAST(ELVART,    3,    1)    ! WRITE CURRENT VALUES
          ELSE
           CALL DATBAST(ELVART,    5,    1)    ! WRITE LAST CONV. VALUES
          ENDIF
          IELEMT=JELEMT
         ENDIF            ! rcoor.ge.-1.0d0....
C
        ENDDO             ! igaust=1,ngault
       ENDIF              ! ndimet.eq.3
C
 2000  CONTINUE
C
      ENDIF               ! nocoit.eq.1
C
 1000 CONTINUE
C
C**** WRITES CONNECTIVITY ARRAY FOR CONTACT ELEMENTS OF NON-COINCIDENT
C     MESH
C
      IF(NOCOIT.GT.0) THEN
       WRITE(LUREST,900)
       DO IELEMT=1,NELEMT
        LGRUPT=MATNOT(IELEMT)
        LTYPET=INT(PROELT( 5,LGRUPT))
        IF(LTYPET.EQ.104) THEN
         NNODNT=INT(PROELT(10,LGRUPT))
         NOCOLT=INT(PROELT(11,LGRUPT))
         IF(NOCOLT.EQ.1) THEN
          WRITE(LUREST,910) IELEMT,MATNOT(IELEMT),
     .                     (LNODST(INODET,IELEMT),INODET=1,NNODNT)
         ENDIF
        ENDIF
       ENDDO
       WRITE(LUREST,902)
  900  FORMAT(//2X,16H CONTACT ELEMENT,5X,5HGROUP,6X,12HNODE NUMBERS)
  902  FORMAT(/)
  910  FORMAT(5X,I7,7X,I5,6X,10I6/10(30X,10I6))
      ENDIF               ! nocoit.eq.1
C
      IF(KERORT.NE.0)
     . CALL RUNENDT('ERROR IN JACOBIAN MATRIX           ')
C
      CALL CPUTIMT(TIME2T)
      CPUSTT=CPUSTT+(TIME2T-TIME1T)
C
      RETURN
      END
