      SUBROUTINE SKYDEST(IDSK2,IFFIX,LNODS,LPNTN,LNUEQ,LPONT,NDOFN,
     .                   NELEM,NEQNS,NLAST,NNODE,NPOIN,NTOTV,NWIDT,
     .                   LEQNS,KRENU)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE COLUMN HEIGHTS FOR THE PROFILE
C
C.... INPUT PARAMETERS
C
c     IFFIX(NDOFN,NPOIN)  -  FIXITY CODES ARRAY
C     LNODS(NNODE,NELEM)  -  ELEMENT CONNECTIVITY
C     LPNTN(NPOIN)        -  RENUMBERING DESTINATIONS
C     NDOFN               -  NO. OF D.O.F PER NODE
C     NELEM               -  NO. OF ELEMENTS
C     NNODE               -  NUMBER OF ELEMENT NODES
C     NPOIN               -  NUMBER OF POINTS
C     NTOTV               -  NUMBER OF VARIABLES
C     NWIDT               -  MAXIMUM BAND-WIDTH ALLOWED
C
C.... OUTPUT PARAMETERS
C
C     LNUEQ(NDOFN,NPOIN)  -  EQUATION NUMBERS
C     LPONT(NEQNS)        -  POSITION OF THE TERM JUST ABOVE THE
C                            DIAGONAL IN THE ARRAYS GSTUP & GSTLO
C     NEQNS               -  NUMBER OF EQUATIONS
C     NLAST               -  NUMBER OF TERMS (TOTAL LENGTH) OF THE
C                            PROFILED ARRAYS GSTUP & GSTLO
C
C.... AUXILIARY WORKING PARAMETERS
C
C     LEQNS(NEVAC)        -  AUXILIARY ARRAY TO STORE THE EQUATION
C                            NUMBERS OF THE CURRENT ELEMENT
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      CHARACTER*80  CAT,CBT,CCT,CDT,CET,CFT,CGT,CHT,CIT,CJT,CKT,CLT,
     .              CMT,CNT,COT,CPT,CQT,CRT,CST,
     .              CA1T,CB1T,CC1T,CD1T,CE1T,CF1T,CG1T,CH1T,CI1T,
     .              C1T,C2T,C3T,C4T,C5T,C6T,C7T,C8T,C9T,C10T
C
      COMMON/NAMETA/CAT,CBT,CCT,CDT,CET,CFT,CGT,CHT,CIT,CJT,CKT,CLT,
     .              CMT,CNT,COT,CPT,CQT,CRT,CST,
     .              CA1T,CB1T,CC1T,CD1T,CE1T,CF1T,CG1T,CH1T,CI1T,
     .              C1T,C2T,C3T,C4T,C5T,C6T,C7T,C8T,C9T,C10T
C
      COMMON/LOGUNT/LUDTST,LUSOLT,LUFROT,LUFRHT,LUDATT,LUPRIT,LUREST,
     .              LUSO2T,LUFR2T,LUPOST,LURSTT,LUBFGT,LUPIPT,LUPANT,
     .              LUGEOT,LUSETT,LUMATT,LUINIT,LULOAT,LUFIXT,LUADVT,
     .              LUACTT,LUFANT,LUSTRT,
     .              LUCU1T,LUCU2T,LUCU3T,LUCU4T,LUCU5T,LUCU6T,LUCU7T,
     .              LUCU8T,LUCU9T,LUC10T
C
      DIMENSION IFFIX(NDOFN,*), LPONT(NTOTV),       LPNTN(*),
     .          LEQNS(*),       LNUEQ(NDOFN,NPOIN), LNODS(NNODE,*)
C
C**** OPEN SOLUTION FILE
C
C     ! ctm: to be revised for general idsk2
C
      IERORT=1
      OPEN(UNIT=LUSOLT,FILE=CBT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     .     ERR=1000)
      IERORT=0
 1000 IF(IERORT.EQ.1) THEN
       WRITE(LUPRIT,901)
       CALL RUNENDT('  ERROR IN OPENING FILES           ')
      ENDIF
  901 FORMAT(' ERROR IN OPENING RESTART FILE 202 (22 LINUX)')
C
      REWIND IDSK2
C
C**** IDENTIFY POSITION AND NUMBER OF EQUATION TO BE SOLVED
C
      IEQNS=0
      IVFIX=0
      DO 5 IPOIN=1,NPOIN
      IF(KRENU.EQ.0) THEN
       IPUNT=IPOIN
      ELSE
       IPUNT=LPNTN(IPOIN)                            ! NODAL RENUMBERING
      ENDIF
C$DIR SCALAR
      DO 5 IDOFN=1,NDOFN
      IF(IFFIX(IDOFN,IPUNT).EQ.0) THEN
        IEQNS=IEQNS+1
        LNUEQ(IDOFN,IPUNT)=IEQNS                     ! NODAL RENUMBERING
      ELSE
        IVFIX=IVFIX+1
        LNUEQ(IDOFN,IPUNT)=-IVFIX
      ENDIF
    5 CONTINUE
C
      NEQNS=IEQNS
C
C**** FIND THE POSITION OF THE DIAGONAL TERMS
C
      NLAST = 0
      DO 10 IEQNS=1,NEQNS
   10 LPONT(IEQNS)=0
C
      IF(NWIDT.LE.0) GO TO 100     ! only diagonal will be assembled
C
      MAXPR = 0
      DO 50 IELEM = 1,NELEM
      MNEQN = 0
      NEVAC = 0
      DO 30 INODE = 1,NNODE
      IPONT = LNODS(INODE,IELEM)
      IF(IPONT.GT.0) THEN
        DO 20 IDOFN = 1,NDOFN
        IEQNS = LNUEQ(IDOFN,IPONT)
        IF(IEQNS.GT.0) THEN
          IF(MNEQN.EQ.0) MNEQN = IEQNS
          MNEQN = MIN0(MNEQN,IEQNS)
          NEVAC = NEVAC + 1
          LEQNS(NEVAC) = IEQNS
        ENDIF
   20   CONTINUE
      ENDIF
   30 CONTINUE
C
      IF(NEVAC.GT.0) THEN
        DO 40 IEVAC = 1,NEVAC
        IEQNS = LEQNS(IEVAC)
        NPONT = MAX0(LPONT(IEQNS),IEQNS-MNEQN)
        LPONT(IEQNS) = MIN0(NPONT,NWIDT)
        MAXPR = MAX0(MAXPR,LPONT(IEQNS))
   40   CONTINUE
      ENDIF
   50 CONTINUE
C
C**** COMPUTE DIAGONAL POINTERS FOR PROFILE
C
      LPONT(1) = 0
      IF(NEQNS.GT.1) THEN
        DO 60 IEQNS = 2,NEQNS
   60   LPONT(IEQNS) = LPONT(IEQNS) + LPONT(IEQNS-1)
        NLAST = LPONT(NEQNS)
      ENDIF
C
  100 IF(NLAST.EQ.0) NLAST=1
      WRITE(IDSK2) NEQNS,NLAST,LNUEQ,LPONT
C
C**** WRITE PROFILE INFORMATION
C
      HEIGM=0
      IF(NEQNS.NE.0) HEIGM=FLOAT(NLAST+NEQNS)/FLOAT(NEQNS)
      MEANH=NINT(HEIGM)
      MAXPR=MAXPR+1
      MAXWD=NWIDT+1
      WRITE(LUREST,900) NTOTV,NEQNS,NLAST,MAXPR,MEANH,MAXWD
      WRITE(LUPRIT,900) NTOTV,NEQNS,NLAST,MAXPR,MEANH,MAXWD
C
      RETURN
  900 FORMAT(//5X,'PROFILE INFORMATION :',/,
     .         5X,'===================  ',/,
     .        15X,'NUMBER OF DOF           =',I8,/
     .        15X,'NUMBER OF EQUATIONS     =',I8,/
     .        15X,'SIZE OF PROFILE         =',I8,' (R*8)'/
     .        15X,'MAXIMUM COLUMN HEIGHT   =',I8,/
     .        15X,'MEAN HALF-WIDTH         =',I8,/
     .        15X,'MAX. HALF-WIDTH ALLOWED =',I8,/)
      END