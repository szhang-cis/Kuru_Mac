      SUBROUTINE INPGEOT(COORDT,ELDATT,INTERT,ITAPET,LNODST,MATNOT,
     .                   PROELT,IBOEMT)
C***********************************************************************
C
C**** THIS ROUTINE READS THE GEOMETRICAL DATA
C
C     Notes: DINT is used to avoid problems in challenge machines
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inpo_omt.f'
C
      DIMENSION ELDATT(*), COORDT(NDIMET,*), LNODST(NNODET,*),
     .          MATNOT(*), PROELT(NPRELT,*)
C
C**** DEFINES MEMORY POINTERS
C
      CALL MEMODTT
C
      IF(IBOEMT.EQ.0) THEN
C
       IERR1T=1
       IF(ITAPET.NE.LUDATT.AND.ITAPET.NE.LUGEOT) GOTO 1000
       IF(ITAPET.EQ.LUGEOT) THEN       ! external .geo file
        OPEN(UNIT=LUGEOT,FILE=CA1T,STATUS='OLD',ERR=1000)
        IERR1T=0
C
 1000   IF(IERR1T.NE.0) THEN
         IF(IERR1T.EQ.1) WRITE(LUREST,901)
         CALL RUNENDT('ERROR IN OPENING FILE         ')
        ENDIF
       ENDIF     ! itapet.eq.lugeot
C
C**** READ THE ELEMENT NODAL CONNECTIONS AND THE PROPERTY NUMBERS
C
       WRITE(LUREST,900)
C
C**** ISTAN DEFINES: =0 standard reading; =1 improved reading
C
       ISTAN=1
C
       NLINET=(NNODET-1)/10+1
       DO 10 IELEMT=1,NELEMT
C
C**** STANDARD READING (ISTAN=0)
C
        IF(ISTAN.EQ.0) THEN
         READ(ITAPET,905) NUMELT,MATNOT(NUMELT),
     .                   (LNODST(INODET,NUMELT),INODET=1,NNODET)
        ELSE
C
C**** IMPROVED READING (ISTAN=1; reading with listent.f)
C
         NPRINT=1
         CALL LISTENT('INPGEOT',NPRINT,ITAPET)
         NUMELT=DINT(PARAMT(1))
         MATNOT(NUMELT)=INT(PARAMT(2))
         DO INODET=1,NNODET
          LNODST(INODET,NUMELT)=INT(PARAMT(2+INODET))
         ENDDO
        ENDIF    ! istan.eq.0
C
C**** WRITE THE ELEMENT NODAL CONECTION AND PROPERTY NUMBERS
C
       WRITE(LUREST,910) NUMELT,MATNOT(NUMELT),
     .                  (LNODST(INODET,NUMELT),INODET=1,NNODET)
   10  CONTINUE
C
C**** ZERO ALL THE NODAL COORDINATES AND THEN READ SOME OF THEM
C        
C$DIR SCALAR
       DO IDIMET=1,NDIMET
        DO IPOINT=1,NPOINT
         COORDT(IDIMET,IPOINT)=0.0D0
        ENDDO
       ENDDO
C
C**** STANDARD READING (ISTAN=0)
C
       IF(ISTAN.EQ.0) THEN  ! to be improved (infinite loop!)
   32   READ(ITAPET,915,ERR=32) 
     .       IPOINT,(COORDT(IDIMET,IPOINT),IDIMET=1,NDIMET)
        IF(IPOINT.NE.NPOINT) GOTO 32
C
C**** IMPROVED READING (ISTAN=1)
C
       ELSE
        DO IPOI1T=1,NPOINT
         CALL LISTENT('INPGEOT',NPRINT,ITAPET)
         IPOINT=INT(PARAMT(1))
         DO IDIMET=1,NDIMET
          COORDT(IDIMET,IPOINT)=PARAMT(IDIMET+1)
         ENDDO
        ENDDO
       ENDIF        ! istan.eq.0
C
C**** LOOP ON ELEMENTS
C
       DO 100 IELEMT=1,NELEMT
       LGRUPT=MATNOT(IELEMT)
C
C**** STORE NUMBER OF NODES PER ELEMENT
C
       NNODLT=0
       DO 40 INODET=1,NNODET
       IPOINT=LNODST(INODET,IELEMT)
  40   IF(IPOINT.NE.0) NNODLT=NNODLT+1
       PROELT(2,LGRUPT)=FLOAT(NNODLT)
C
C**** INTERPOLATE MID-SIDE NODES IF NECESSARY
C
       IF(INTERT.EQ.1) THEN
        IF(NDIMET.EQ.1) THEN
         CALL NODES1(COORDT,LNODST(1,IELEMT),NDIMET,IELEMT,
     .               NNODLT)
        ELSE IF(NDIMET.EQ.2.OR.KPROBT.EQ.2) THEN
         CALL NODES2(COORDT,LNODST(1,IELEMT),NDIMET,IELEMT,
     .               NNODLT)
        ELSE IF(NDIMET.EQ.3) THEN
         CALL NODES3(COORDT,LNODST(1,IELEMT),NDIMET,IELEMT,
     .               NNODLT)
        ENDIF
       ENDIF
C
  100  CONTINUE                 ! END LOOP ON ELEMENTS
C
       IF(NMEMO1.EQ.0) THEN     ! coordinates in an elemental array
C
C**** STORE ELEMENT COORDINATES IN DATA-BASE
C
        DO IELEMT=1,NELEMT
         LGRUPT=MATNOT(IELEMT)
         NNODLT=INT(PROELT(2,LGRUPT))
          CALL GATHER(COORDT,NDIMET,NPOINT,ELDATT(IDATAT(1)),NDIMET,
     .                NNODLT,LNODST(1,IELEMT))
          CALL DATBAST(ELDATT,    1,    1)
        ENDDO
       ENDIF
C
C**** PRINTOUT NODAL COORDINATES 
C
       IF(NDIMET.EQ.2) THEN 
        WRITE(LUREST,950)
        WRITE(LUREST,971)
     .   (IPOINT,(COORDT(IDIMET,IPOINT),IDIMET=1,NDIMET),
     .    IPOINT=1,NPOINT)
       ELSE IF(NDIMET.EQ.3) THEN
        IF(KPROBT.NE.2) THEN                 ! not shells
         WRITE(LUREST,960)
         WRITE(LUREST,970)(IPOINT,(COORDT(IDIMET,IPOINT),IDIMET=1,
     .                     NDIMET),IPOINT=1,NPOINT)
        ELSE
          WRITE(LUREST,980)
          WRITE(LUREST,970)(IPOINT,(COORDT(IDIMET,IPOINT),
     .                      IDIMET=1,NDIMET),
     .                     (COORDT(IDIMET,IPOINT+NPOINT),
     .                      IDIMET=1,NDIMET),IPOINT=1,NPOINT)
        ENDIF
       ENDIF
C
       WRITE(LUREST,902)
C
C**** LOOK FOR 'END_GEOMETRY' CARD
C
       NPRINT=0
       JTAPET=LUDATT
       CALL LISTENT('INPGEOT',NPRINT,JTAPET)
       IF(WORDST(1).NE.'END_G')
     .  CALL RUNENDT('INPGEOT:END_GEOMETRY CARD NOT FOUND')
       WRITE(LUREST,902)
C
      ELSE              ! iboemt=1
C
C**** CHECK THE INPUT DATA
C
       CALL INPCEKT(COORDT,LNODST,MATNOT,NDIMET,NELEMT,NMATST,NNODET,
     .              NPOINT,PROELT,NPRELT,NGRUPT,NOCOIT,LUREST)
C
      ENDIF             ! iboemt.eq.0
C
      RETURN
C
  900 FORMAT(//2X,8H ELEMENT,5X,5HGROUP,6X,12HNODE NUMBERS)
  901 FORMAT(' ERROR IN OPENING GEOMETRY INPUT FILE 240 (50 LINUX)')
  902 FORMAT(/)
  905 FORMAT(12I5/10(10X,10I5))
  910 FORMAT(5X,I7,5X,I5,6X,10I6/10(26X,10I6))
  915 FORMAT(I5,6F15.6)
  950 FORMAT(//3(5X,5H NODE,2X,8H X-COORD,2X,8H Y-COORD))
  960 FORMAT(//3(5X,5H NODE,2X,8H X-COORD,2X,8H Y-COORD,2X,8H Z-COORD))
  970 FORMAT(3(5X,I6,3F10.3))
  971 FORMAT(3(5X,I6,2F10.3))
  980 FORMAT(//(5X,5H NODE,2X,8H X-COORD,2X,8H Y-COORD,2X,8H Z-COORD,
     .                      2X,8H  ALFA-1,2X,8H  ALFA-2,2X,8H  ALFA-3))
c  990 FORMAT(5X,I5,6F10.3)
      END