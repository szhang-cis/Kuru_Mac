      SUBROUTINE SETMTX(ELDAT,ELPRE,ELVAR,ELMAT,LNODS,MATNO,PROEL,PROPS,
     .                  COORD,VNORM,DISTO,SINI2,SINI3,WORK1,IPRIX)
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
      INCLUDE 'addi_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
      INCLUDE 'inpo_om.f'
C
      COMMON/JACOBSA/IEROR,KEROR
C
      DIMENSION MATNO(NELEM),       LNODS(NNODE,NELEM),
     .          PROEL(NPREL,NGRUP), PROPS(NPROP,NMATS),
     .          ELDAT(NDATA),       ELPRE(NPREV),
     .          ELVAR(NSTAT),       ELMAT(NMATX),
     .          COORD(NDIME,NPOIN), VNORM(*),
     .          DISTO(*),           WORK1(*)
      DIMENSION SINI2(*),           SINI3(*)
C
      CALL CPUTIM(TIME1)
C
C**** LOOP ON ELEMENTS
C
C     Note: KELEM is used instead of IELEM as a loop index in order to
C           avoid compiling problems when NOCOI>0 is considered since
C           in this case a JELEM do loop is also performed and IELEM
C           is reserved to be used in datbast.f for these two loops
C           (IELEM is transferred by auxl_om.f).
C
      KEROR=0
      DO 1000 KELEM=1,NELEM
      IELEM=KELEM
      LGRUP=MATNO(KELEM)
      LMATS=INT(PROEL(1,LGRUP))
      NNODL=INT(PROEL(2,LGRUP))
      LTYPE=INT(PROEL(5,LGRUP))
C
      IEROR=0
C
      IF(NMEMO1M.EQ.0) THEN          ! coordinates in an elemental array
C
C**** READ ELDAT ( coordinates ) FROM DATA-BASE
C
       CALL DATBAS(ELDAT,    1,    2)
      ELSE                           ! coordinates in a global array
       CALL GATHER(COORD,NDIME,NPOIN,WORK1(ISETM(10)),NDIME,NNODL,
     .             LNODS(1,KELEM))
      ENDIF
C
      IF(NOCOI.GT.0) THEN
       IF(LTYPE.EQ.32) THEN
        NOCOL=INT(PROEL(22,LGRUP))
        IF(NOCOL.EQ.1) THEN
         IF(LARGC.NE.0) THEN
          IF(NMEMO5M.EQ.0) THEN
C
C**** SCALAR GATHER OPERATIONS ( DISTO ---> ELDIS(ELVAR) )
C
           CALL GATHER(DISTO,NDOFC,NPOIN,ELVAR,           NDOFC,NNODL,
     .                 LNODS(1,KELEM))
         ELSE
C
C**** SCALAR GATHER OPERATIONS ( DISTO ---> ELDIS )
C
           CALL GATHER(DISTO,NDOFC,NPOIN,WORK1(ISETM(17)),NDOFC,NNODL,
     .                 LNODS(1,KELEM))
          ENDIF                  ! nmemo5m.eq.0
         ENDIF                   ! largc.ne.0
        ENDIF                    ! nocol.eq.1
       ENDIF                     ! ltype.eq.32
      ENDIF                      ! nocoi.gt.0
C
      INITX=0
      IF(INITV.EQ.1.AND.IPRCO.EQ.0) THEN    ! non-standard initial cond.
       IF(LTYPE.EQ.30) THEN
        IF(NPRE2.GT.0) THEN
C
C**** SCALAR GATHER OPERATIONS ( SINI2 ---> WORK1(ISETM(32))    )
C
         CALL GATHER(SINI2,NPRE2,NPOIN,WORK1(ISETM(32)),NPRE2,
     .               NNODL,LNODS(1,IELEM))
C
C**** READ CURRENT PRESCRIBED VARIABLES FROM DATA BASE
C
         CALL DATBAS(ELPRE,    2,    2)
        ENDIF
        IF(NPRE3.GT.0) THEN
C
C**** SCALAR GATHER OPERATIONS ( SINI3 ---> WORK1(ISETM(33))    )
C
         CALL GATHER(SINI3,NPRE3,NPOIN,WORK1(ISETM(33)),NPRE3,
     .               NNODL,LNODS(1,IELEM))
C
C**** READ CURRENT STATE VARIABLES FROM DATA BASE
C
         CALL DATBAS(ELVAR,    3,    2)
         INITX=1
        ENDIF
       ENDIF
      ENDIF
C
C**** READ CURRENT STATE VARIABLES FROM DATA BASE
C
      IF(INITP.EQ.1.AND.IPRCO.EQ.0.AND.INITX.EQ.0) THEN
       CALL DATBAS(ELVAR,    3,    2)
      ENDIF
C
C**** SET UP ARRAY ELDAT & WRITE IT TO DATA-BASE
C
      CALL ELMLIB(LNODS(1,KELEM),PROEL(1,LGRUP),PROPS(1,LMATS),
     .            WORK1,WORK1,WORK1,WORK1,WORK1,
     .            ELDAT,ELPRE,ELVAR,ELMAT,WORK1,    1)
C
      IF(IEROR.NE.0) GO TO 1000
C
C**** WRITE ELDAT TO DATA-BASE
C
      IF(NMEMO1M.EQ.0.OR.NMEMO2M.EQ.0)
     . CALL DATBAS(ELDAT,    1,    1)
C
      IF(INITV.EQ.1.AND.IPRCO.EQ.0) THEN    ! non-standard initial cond.
       IF(LTYPE.EQ.30) THEN
C
C**** WRITE CURRENT PRESCRIBED & STATE VARIABLES TO DATA BASE
C     (see adddat.f)
C
        IF(NPRE2.GT.0) CALL DATBAS(ELPRE,    2,    1)
        IF(NPRE3.GT.0) CALL DATBAS(ELVAR,    3,    1)
       ENDIF
      ENDIF
C
C**** WRITE CURRENT STATE VARIABLES TO DATA BASE
C
      IF(INITP.EQ.1.AND.IPRCO.EQ.0.AND.INITX.EQ.0) THEN
       CALL DATBAS(ELVAR,    3,    1)
      ENDIF
C
C**** SCALAR SCATTER OPERATION ( WORK1 --> VNORM )
C
      CALL SCATER(WORK1(ISETM(9)),NDIME,NNODL,VNORM,NDOFC,NPOIN,
     .            LNODS(1,KELEM))
C
C**** DEALS WITH CONTACT ELEMENTS FOR NON-COINCIDENT MESH
C
      IF(NOCOI.GT.0) THEN
       IF(LTYPE.NE.32) GO TO 1000            ! skip non-contact elements
       IF(NOCOL.EQ.0) GO TO 1000             ! skip coincident cont. el.
C
       NGAUL=INT(PROEL(4,LGRUP))
C
       ISETC=INT(PROEL(19,LGRUP))            ! contact set
       IMAES=INT(PROEL(20,LGRUP))            ! contact index
C
C**** RECOVERS COORDINATES & NORMAL OF MASTER CONTACT ELEMENT (KELEM)
C     AT GAUSS POINTS
C
       IF(NMEMO1M.EQ.0) THEN
        IAUXX=-1
        DO IGAUS=1,NGAUL
         DO IDIME=1,NDIME
          IAUXX=IAUXX+1
          WORK1(ISETM(13)+IAUXX)=ELDAT(IDATA( 4)+IAUXX)          ! GPCOD
          WORK1(ISETM(23)+IAUXX)=ELDAT(IDATA(11)+IAUXX)          ! VNORL
         ENDDO
        ENDDO
       ELSE
        IAUXX=-1
        DO IGAUS=1,NGAUL
         DO IDIME=1,NDIME
          IAUXX=IAUXX+1
          WORK1(ISETM(21)+IAUXX)=WORK1(ISETM(13)+IAUXX)          ! GPCOD
          WORK1(ISETM(23)+IAUXX)=WORK1(ISETM( 9)+IAUXX)          ! VNORL
         ENDDO
        ENDDO
       ENDIF
C
       DO IGAUS=1,NGAUL               ! proper initialization for solver
        DO INODL=1,NNODL
         INODX=NNODL*IGAUS+INODL
         LNODS(NNODL*IGAUS+INODL,KELEM)=LNODS(INODL,KELEM)
        ENDDO
       ENDDO
C
       DO IGAUS=1,NGAUL
        WORK1(ISETM(31)+IGAUS-1)=-1.0D+10 ! GAPMI: min. distance along n
       ENDDO
C
C**** LOOP ON ELEMENTS
C
C     Note: this loop is only necessary to compute nodal coordinates of
C           the slave elements
C
       DO 2000 JELEM=1,NELEM
       IF(KELEM.EQ.JELEM) GO TO 2000         ! skip equal elements
       IELEM=JELEM
       JGRUP=MATNO(JELEM)
       JMATS=INT(PROEL(1,JGRUP))
       JNODL=INT(PROEL(2,JGRUP))
       JTYPE=INT(PROEL(5,JGRUP))
       IF(JTYPE.NE.32) GO TO 2000            ! skip non-contact elements
       JOCOL=INT(PROEL(22,JGRUP))
       IF(JOCOL.EQ.0) GO TO 2000             ! skip coincident cont. el.
C
       JSETC=INT(PROEL(19,JGRUP))            ! contact set
       JMAES=INT(PROEL(20,JGRUP))            ! contact index
C
C**** CONTACT CRITERIA
C
       IAUXX=0
       IF((IMAES.EQ.1.AND.JMAES.EQ.1).OR.    ! no self-contact (nsc)
     .    (IMAES.EQ.1.AND.JMAES.EQ.3).OR.
     .    (IMAES.EQ.3.AND.JMAES.EQ.1)) THEN
        IF(ISETC.EQ.JSETC) THEN
         IF(LGRUP.NE.JGRUP) IAUXX=1
        ENDIF
       ENDIF
       IF((IMAES.EQ.2.AND.JMAES.EQ.2).OR.    ! self-contact (sc)
     .    (IMAES.EQ.2.AND.JMAES.EQ.3).OR.
     .    (IMAES.EQ.3.AND.JMAES.EQ.2)) THEN
        IF(ISETC.EQ.JSETC) THEN
         IF(LGRUP.EQ.JGRUP) IAUXX=1
        ENDIF
       ENDIF
       IF(IMAES.EQ.3.AND.JMAES.EQ.3) THEN    ! both (nsc & sc)
        IF(ISETC.EQ.JSETC) IAUXX=1
       ENDIF
       IF(IAUXX.EQ.0) GO TO 2000          ! skip due to contact criteria
C
       IF(NMEMO1M.EQ.0) THEN    ! coordinates in an elemental array
C
C**** READ ELDAT ( coordinates ) FROM DATA-BASE
C
        CALL DATBAS(ELDAT,    1,    2)
C
       ELSE                     ! coordinates in a global array
        CALL GATHER(COORD,NDIME,NPOIN,WORK1(ISETM(10)),NDIME,
     .              NNODL,LNODS(1,JELEM))
       ENDIF
C
       IF(LARGC.NE.0) THEN
        IF(NMEMO5M.EQ.0) THEN
C
C**** SCALAR GATHER OPERATIONS ( DISTO ---> ELDIS(ELVAR) )
C
         CALL GATHER(DISTO,NDOFC,NPOIN,ELVAR,           NDOFC,NNODL,
     .               LNODS(1,JELEM))
        ELSE
C
C**** SCALAR GATHER OPERATIONS ( DISTO ---> ELDIS )
C
         CALL GATHER(DISTO,NDOFC,NPOIN,WORK1(ISETM(17)),NDOFC,NNODL,
     .               LNODS(1,JELEM))
        ENDIF                   ! nmemo5m.eq.0
       ENDIF                    ! largc.ne.0
C
C**** SET UP ARRAY ELDAT & WRITE IT TO DATA BASE
C
       CALL ELMLIB(LNODS(1,JELEM),PROEL(1,JGRUP),PROPS(1,JMATS),
     .             WORK1,WORK1,WORK1,WORK1,WORK1,
     .             ELDAT,ELPRE,ELVAR,ELMAT,WORK1,    1)
C
C**** RECOVERS COORDINATES OF SLAVE CONTACT ELEMENT (JELEM) AT NODES
C
       IF(NMEMO1M.EQ.0) THEN
        IAUXX=-1
        DO IGAUS=1,NGAUL
         DO IDIME=1,NDIME
          IAUXX=IAUXX+1
          WORK1(ISETM(22)+IAUXX)=ELDAT(IDATA( 1)+IAUXX)          ! ELCOD
         ENDDO
        ENDDO
       ELSE
        IAUXX=-1
        DO INODL=1,NNODL
         DO IDIME=1,NDIME
          IAUXX=IAUXX+1
          WORK1(ISETM(22)+IAUXX)=WORK1(ISETM(10)+IAUXX)          ! ELCOD
         ENDDO
        ENDDO
       ENDIF
C
C**** TASKS (ACCORDING TO NDIME):
C
C     1) COMPUTES COORDINATES OF SLAVE CONTACT POINT CORRESPONDING TO
C        THE GAUSS POINT OF THE MASTER ELEMENT (THESE TWO POINTS ARE
C        THOSE PRESENTING THE MINIMUN DISTANCE BETWEEN THEM)
C
C     2) REDEFINES CONNECTIVITY ARRAY OF THE MASTER ELEMENT BY ADDING
C        NODES OF THE SLAVE ELEMEMT
C
C     3) COMPUTES NATURAL COORDINATES OF CONTACT POINT & ASSIGNS THEM TO
C        ELVAR (STRSG)
C
C
C        Notes:
C
C        Index IAUXZ of STRSG(NDIME-1): natural coordinates of
C        contact point
C
       IF(NDIME.EQ.2) THEN
        X1A=WORK1(ISETM(22)  )
        Y1A=WORK1(ISETM(22)+1)
        X1B=WORK1(ISETM(22)+2)
        Y1B=WORK1(ISETM(22)+3)
        DO IGAUS=1,NGAUL
         X1= WORK1(ISETM(21)+(IGAUS-1)*NDIME  )
         Y1= WORK1(ISETM(21)+(IGAUS-1)*NDIME+1)
         PX1=WORK1(ISETM(23)+(IGAUS-1)*NDIME  )
         PY1=WORK1(ISETM(23)+(IGAUS-1)*NDIME+1)
C
         CALL CPCALC(X1A,Y1A,Z1A,X1B,Y1B,Z1B,X1C,Y1C,Z1C,X1D,Y1D,Z1D,
     .               X1, Y1, Z1, PX1,PY1,PZ1,
     .               NDIME,NNODL,RCOOR,SCOOR,DISTT,DISMI,IAUXE,
     .               RCINF,SCINF,RCSUP)
         IF(IAUXE.EQ.1)
     .    CALL RUNEND('ERROR: NATURAL COORD. NOT FOUND IN cpcalc.f')
C
         GAPMI=WORK1(ISETM(31)+IGAUS-1)
         IF((RCOOR.GE.-1.0D0.AND.RCOOR.LE.1.0D0).AND.
     .                       (DISTT.LT.GAPPC).AND.(DISMI.GT.GAPMI)) THEN
          WORK1(ISETM(31)+IGAUS-1)=DISMI
          IAUXY=(IGAUS-1)*NNODL
          LNODS(NNODL+IAUXY+1,KELEM)=LNODS(1,JELEM)
          LNODS(NNODL+IAUXY+2,KELEM)=LNODS(2,JELEM)
C
          IELEM=KELEM
          IF(IPRCO.EQ.0) THEN
           CALL DATBAS(ELVAR,    3,    2)      ! READ CURRENT VALUES
          ELSE
           CALL DATBAS(ELVAR,    5,    2)      ! READ LAST CONV. VALUES
          ENDIF
C
          IAUXZ=(IGAUS-1)*NSTR1
          ELVAR(ISTAT(4)+IAUXZ)=RCOOR
C
          IF(IPRCO.EQ.0) THEN
           CALL DATBAS(ELVAR,    3,    1)      ! WRITE CURRENT VALUES
          ELSE
           CALL DATBAS(ELVAR,    5,    1)      ! WRITE LAST CONV. VALUES
          ENDIF
          IELEM=JELEM
         ENDIF            ! rcoor.ge.-1.0d0....
C
        ENDDO             ! igaus=1,ngaul
       ENDIF              ! ndime.eq.2
C
       IF(NDIME.EQ.3) THEN
        X1A=WORK1(ISETM(22)  )
        Y1A=WORK1(ISETM(22)+1)
        Z1A=WORK1(ISETM(22)+2)
        X1B=WORK1(ISETM(22)+3)
        Y1B=WORK1(ISETM(22)+4)
        Z1B=WORK1(ISETM(22)+5)
        X1C=WORK1(ISETM(22)+6)
        Y1C=WORK1(ISETM(22)+7)
        Z1C=WORK1(ISETM(22)+8)
        IF(NNODL.GT.3) THEN
         X1D=WORK1(ISETM(22)+ 9)
         Y1D=WORK1(ISETM(22)+10)
         Z1D=WORK1(ISETM(22)+11)
        ENDIF
        DO IGAUS=1,NGAUL
         X1= WORK1(ISETM(21)+(IGAUS-1)*NDIME  )
         Y1= WORK1(ISETM(21)+(IGAUS-1)*NDIME+1)
         Z1= WORK1(ISETM(21)+(IGAUS-1)*NDIME+2)
         PX1=WORK1(ISETM(23)+(IGAUS-1)*NDIME  )
         PY1=WORK1(ISETM(23)+(IGAUS-1)*NDIME+1)
         PZ1=WORK1(ISETM(23)+(IGAUS-1)*NDIME+2)
C
         CALL CPCALC(X1A,Y1A,Z1A,X1B,Y1B,Z1B,X1C,Y1C,Z1C,X1D,Y1D,Z1D,
     .               X1, Y1, Z1, PX1,PY1,PZ1,
     .               NDIME,NNODL,RCOOR,SCOOR,DISTT,DISMI,IAUXE,
     .               RCINF,SCINF,RCSUP)
         IF(IAUXE.EQ.1)
     .    CALL RUNEND('ERROR: NATURAL COORD. NOT FOUND IN cpcalc.f')
C
         GAPMI=WORK1(ISETM(31)+IGAUS-1)
         IF((RCOOR.GE.RCINF.AND.RCOOR.LE.1.0D0).AND.
     .      (SCOOR.GE.SCINF.AND.SCOOR.LE.(1.0D0-RCSUP)).AND.
     .                       (DISTT.LT.GAPPC).AND.(DISMI.GT.GAPMI)) THEN
          WORK1(ISETM(31)+IGAUS-1)=DISMI
          IAUXY=(IGAUS-1)*NNODL
          LNODS(NNODL+IAUXY+1,KELEM)=LNODS(1,JELEM)
          LNODS(NNODL+IAUXY+2,KELEM)=LNODS(2,JELEM)
          LNODS(NNODL+IAUXY+3,KELEM)=LNODS(3,JELEM)
          IF(NNODL.GT.3)
     .    LNODS(NNODL+IAUXY+4,KELEM)=LNODS(4,JELEM)
C
          IELEM=KELEM
          IF(IPRCO.EQ.0) THEN
           CALL DATBAS(ELVAR,    3,    2)      ! READ CURRENT VALUES
          ELSE
           CALL DATBAS(ELVAR,    5,    2)      ! READ LAST CONV. VALUES
          ENDIF
C
          IAUXZ=(IGAUS-1)*NSTR1
          ELVAR(ISTAT(4)+IAUXZ  )=RCOOR
          ELVAR(ISTAT(4)+IAUXZ+1)=SCOOR
C
          IF(IPRCO.EQ.0) THEN
           CALL DATBAS(ELVAR,    3,    1)      ! WRITE CURRENT VALUES
          ELSE
           CALL DATBAS(ELVAR,    5,    1)      ! WRITE LAST CONV. VALUES
          ENDIF
          IELEM=JELEM
         ENDIF            ! rcoor.ge.-1.0d0....
C
        ENDDO             ! igaus=1,ngaul
       ENDIF              ! ndime.eq.3
C
 2000  CONTINUE
C
      ENDIF               ! nocoi.gt.0
C
 1000 CONTINUE
C
C**** WRITES CONNECTIVITY ARRAY FOR CONTACT ELEMENTS OF NON-COINCIDENT
C     MESH
C
      IF(NOCOI.GT.0.AND.IPRIX.EQ.1) THEN
       WRITE(LURES,900)
       DO IELEM=1,NELEM
        LGRUP=MATNO(IELEM)
        LTYPE=INT(PROEL( 5,LGRUP))
        IF(LTYPE.EQ.32) THEN
         NNODN=INT(PROEL(21,LGRUP))
         NOCOL=INT(PROEL(22,LGRUP))
         IF(NOCOL.EQ.1) THEN
          WRITE(LURES,910) IELEM,MATNO(IELEM),
     .                    (LNODS(INODE,IELEM),INODE=1,NNODN)
         ENDIF
        ENDIF
       ENDDO
       WRITE(LURES,902)
  900  FORMAT(//2X,16H CONTACT ELEMENT,5X,5HGROUP,6X,12HNODE NUMBERS)
  902  FORMAT(/)
  910  FORMAT(5X,I7,7X,I5,6X,10I6/10(30X,10I6))
      ENDIF               ! nocoi.gt.0.and.iprix.eq.1
C
      IF(KEROR.NE.0)
     . CALL RUNEND('ERROR IN JACOBIAN MATRIX           ')
C
      CALL CPUTIM(TIME2)
      CPUST=CPUST+(TIME2-TIME1)
C
      RETURN
      END
