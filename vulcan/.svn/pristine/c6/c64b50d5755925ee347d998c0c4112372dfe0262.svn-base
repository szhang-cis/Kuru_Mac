      SUBROUTINE INPCEKS(COORD,LNODS,MATNO,NDIME,NELEM,NMATS,NNODE,
     .                   NPOIN,proel,nprel,ngrup)
C***********************************************************************
C
C**** THIS ROUTINE CHECKS THE GEOMETRICAL INPUT DATA
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
c     COMMON/LOGUNT/LUDTST,LUSOLT,LUFROT,LUFRHT,LUDATT,LUPRIT,LUREST,
c    .              LUSO2T,LUFR2T,LUPOST,LURSTT,LUBFGT,LUPIPT,LUPANT,
c    .              LUGEOT,LUSETT,LUMATT,LUINIT,LULOAT,LUFIXT,LUADVT,
c    .              LUACTS,LUFANS,
c    .              LUCU1T,LUCU2T,LUCU3T,LUCU4T,LUCU5T,LUCU6T,LUCU7T,
c    .              LUCU8T,LUCU9T,LUC10T
      COMMON/LOGUNS/LUSOLS,LUFRHS,LUDATS,LUPRIS,LURESS,LUPOSS,
     .              LUCU1S
C
      DIMENSION COORD(NDIME,*), LNODS(NNODE,*),
     .          MATNO(*)
      DIMENSION NEROR(24),      proel(nprel,ngrup)
C
C**** INITIALISED THE LIST OF ERRORS
C
      DO 5 IEROR=13,24
    5 NEROR(IEROR)=0
C
C**** CHECK AGAINST TWO IDENTICAL NONZERO NODAL COORDINATES
C
      DO 40 IPOIN=2,NPOIN
C$DIR SCALAR
      DO 30 JPOIN=1,IPOIN-1
C$DIR SCALAR
      DO 20 IDIME=1,NDIME
   20 IF(COORD(IDIME,IPOIN).NE.COORD(IDIME,JPOIN)) GOTO 30
C
C**** FORM TO AVOID CHECK OF COINCIDENT COORDINATES FOR THERMAL
C     CONTACT ELEMENTS
C
      jjj=0
      do ielem=1,nelem
       igrup=matno(ielem)
       itype=int(proel(5,igrup))
       NNODL=INT(PROEL(2,IGRUP))
       if(itype.eq.104) then
C
        if(ndime.eq.1) then
         if(nnodl.eq.2) then
          iii=0
          do inode=1,nnode
           if(lnods(inode,ielem).eq.ipoin) iii=iii+1
           if(lnods(inode,ielem).eq.jpoin) iii=iii+1
          enddo
          if(iii.eq.2) go to 30     ! same element
          if(iii.eq.1) then         ! different element
           jjj=jjj+1
           if(jjj.eq.2) go to 30
          endif
         else
          CALL RUNENDS('INPCEK:ERROR DETECTED WITH NNODL.NE.2 ')
         endif
        endif              ! ndime.eq.1
C
        if(ndime.eq.2) then
         if(nnodl.eq.4) then
          iii=0
          inodi=0
          inodj=0
          do inode=1,nnode
           if(lnods(inode,ielem).eq.ipoin) then
            iii=iii+1
            inodi=inode
           endif
           if(lnods(inode,ielem).eq.jpoin) then
            iii=iii+1
            inodj=inode
           endif
          enddo
          if(iii.eq.2) then         ! same element
           if(inodi.eq.0.or.inodj.eq.0)
     .      call runends('inodi-inodj = 0       ')
           if(inodi.eq.1.and.inodj.eq.4) go to 30
           if(inodi.eq.4.and.inodj.eq.1) go to 30
           if(inodi.eq.2.and.inodj.eq.3) go to 30
           if(inodi.eq.3.and.inodj.eq.2) go to 30
          endif
          if(iii.eq.1) then         ! different element
           jjj=jjj+1
           if(jjj.eq.2) go to 30
          endif
         else
          CALL RUNENDS('INPCEK:ERROR DETECTED WITH NNODL.NE.4 ')
         endif
        endif              ! ndime.eq.2
C
        if(ndime.eq.3) then
         if(nnodl.eq.6.or.nnodl.eq.8) then
          iii=0
          inodi=0
          inodj=0
          do inode=1,nnode
           if(lnods(inode,ielem).eq.ipoin) then
            iii=iii+1
            inodi=inode
           endif
           if(lnods(inode,ielem).eq.jpoin) then
            iii=iii+1
            inodj=inode
           endif
          enddo
          if(iii.eq.2) then         ! same element
           if(inodi.eq.0.or.inodj.eq.0)
     .      call runends('inodi-inodj = 0       ')
           if(nnodl.eq.6) then
            if(inodi.eq.1.and.inodj.eq.4) go to 30
            if(inodi.eq.4.and.inodj.eq.1) go to 30
            if(inodi.eq.2.and.inodj.eq.5) go to 30
            if(inodi.eq.5.and.inodj.eq.2) go to 30
            if(inodi.eq.3.and.inodj.eq.6) go to 30
            if(inodi.eq.6.and.inodj.eq.3) go to 30
           endif
           if(nnodl.eq.8) then
            if(inodi.eq.1.and.inodj.eq.5) go to 30
            if(inodi.eq.5.and.inodj.eq.1) go to 30
            if(inodi.eq.2.and.inodj.eq.6) go to 30
            if(inodi.eq.6.and.inodj.eq.2) go to 30
            if(inodi.eq.3.and.inodj.eq.7) go to 30
            if(inodi.eq.7.and.inodj.eq.3) go to 30
            if(inodi.eq.4.and.inodj.eq.8) go to 30
            if(inodi.eq.8.and.inodj.eq.4) go to 30
           endif
          endif
          if(iii.eq.1) then         ! different element
           jjj=jjj+1
           if(jjj.eq.2) go to 30
          endif
         else
          CALL RUNENDS('INPCEK:ERROR DETECTED WITH NNODL.NE.6-8 ')
         endif
        endif              ! ndime.eq.3
C
       endif
      enddo
C
      NEROR(13)=NEROR(13)+1
      WRITE(LURESS,900)IPOIN,JPOIN
   30 CONTINUE
   40 CONTINUE
C
C**** COMPATIBLE INPUT WITH THE REST OF ELEMENTS
C
C     INPUT 4-----3         OUTPUT 3-----4
C           .     .                .     .
C           1-----2                1-----2
C
      if(ndime.eq.2) then
       do ielem=1,nelem
        igrup=matno(ielem)
        itype=int(proel(5,igrup))
        NNODL=INT(PROEL(2,IGRUP))
        if(itype.eq.104) then
         if(nnodl.eq.4) then
          lauxi=lnods(3,ielem)
          lnods(3,ielem)=lnods(4,ielem)
          lnods(4,ielem)=lauxi
         else
          CALL RUNENDS('INPCEK:ERROR DETECTED WITH NNODL.NE.4 ')
         endif
        endif
       enddo
      endif
C
C**** CHECK THE LIST OF ELEMENT PROPERTY NUMBERS
C
C$DIR SCALAR
      DO 50 IELEM=1,NELEM
      IF((MATNO(IELEM).LE.0).OR.(MATNO(IELEM).GT.NMATS)) THEN
       NEROR(14)=NEROR(14)+1
       WRITE(LURESS,905)IELEM,MATNO(IELEM)
      ENDIF
   50 CONTINUE
C
C**** CHECK FOR IMPOSSIBLE NODE NUMBERS
C
      DO 60 IELEM=1,NELEM
C$DIR SCALAR
      DO 60 INODE=1,NNODE
      IF((LNODS(INODE,IELEM).LT.0).OR.
     .   (LNODS(INODE,IELEM).GT.NPOIN)) THEN
       NEROR(16)=NEROR(16)+1
       WRITE(LURESS,910)IELEM,INODE,LNODS(INODE,IELEM),NPOIN
      ENDIF
   60 CONTINUE
C
C**** CHECK FOR ANY REPETITION OF A NODE NUMBER WITHIN AN ELEMENT
C
      NPOIT=NPOIN
      NPOIT=0                         ! >>>>> TOO MUCH WORK 
      DO 110 IPOIN=1,NPOIT
C
C**** SEEK FIRST,LAST AND INTERMEDIATE APPEARANCES OF NODE IPOIN & 
C     CALCULATE INCREASE OR DECREASE IN FRONTWIDTH AT EACH ELEMENT
C     STAGE
C
      KSTAR=0
      DO 120 IELEM=1,NELEM
      KZERO=0
      DO 130 INODE=1,NNODE
      IF(LNODS(INODE,IELEM).NE.IPOIN)GOTO 130
      KZERO=KZERO+1
C
      IF(KZERO.GT.1)THEN
       NEROR(17)=NEROR(17)+1
       WRITE(LURESS,915)IPOIN,IELEM
      ENDIF
C
      IF(KSTAR.EQ.0)KSTAR=IELEM
C
      KLAST=IELEM
      NLAST=INODE
  130 CONTINUE
  120 CONTINUE
C
C**** CHECK IF THIS IS AN UNUSED NODE & IF IT HAS NON-ZERO COORDINATES
C
      IF(KSTAR.EQ.0)THEN
       WRITE(LURESS,920) IPOIN
       NEROR(18)=NEROR(18)+1
C
       SIGMA=0.0
C$DIR SCALAR
       DO 140 IDIME=1,NDIME
  140  SIGMA=SIGMA+COORD(IDIME,IPOIN)*COORD(IDIME,IPOIN)
       IF(SIGMA.NE.0.)THEN
        NEROR(19)=NEROR(19)+1
        WRITE(LURESS,925)
       ENDIF
C
      ENDIF
C
  110 CONTINUE
C
C**** CHECKS DIFFERENT SIGN IN THE R COORDINATE FOR AXISYMMETRY
C
      IPOSI=0
      INEGA=0
      DO IELEM=1,NELEM
       LGRUP=MATNO(IELEM)
       NTYPE=INT(PROEL(6,LGRUP))
       NNODL=INT(PROEL(2,LGRUP))
       IF(NTYPE.EQ.3) THEN
        DO INODL=1,NNODL
         IPOIN=LNODS(INODL,IELEM)
         IF(COORD(1,IPOIN).GE.0.0) THEN
          IPOSI=IPOSI+1
         ELSE
          INEGA=INEGA+1
         ENDIF
        ENDDO
       ENDIF
      ENDDO
      IF(IPOSI.GT.0.AND.INEGA.GT.0)
     . CALL RUNENDS('ERROR: WRONG R COORDINATES')
C
C**** VERIFIED IF ANY INPUT ERRORS HAVE BEEN DETECTED
C
      KEROR=0
      DO 200 IEROR=13,24
  200 IF(NEROR(IEROR).NE.0) KEROR=1
      IF(KEROR.EQ.1)  
     . CALL RUNENDS('ERROR DETECTED IN INPCEK           ')
C
      RETURN
  900 FORMAT(1X,
     .'IDENTICAL COORDINATES HAS BEEN FOUND FOR NODAL POINTS N0.',2I5)
  905 FORMAT(1X,'MATERIAL NO. OUT OF RANGE 1-NMATS. MATNO(',I5,')=',I5) 
  910 FORMAT(11X,'IMPOSSIBLE NODE NUMBER LNODS(',I5,',',I5,')=',I5,/,
     .'SINCE IS EITHER LESS THAN ZERO OR GREATER THAN THE MAXIMUM',/,
     .'SPECIFIED NUMBER OF POINTS.  NPOIN=',I5)
  915 FORMAT(/,
     .' POINT NO.',I5,' APPEARS MORE THAN ONES IN THE LIST OF',/,
     .' NODAL CONNECTIONS OF ELEMENT NO.',I5)
  920 FORMAT(/,15H CHECK WHY NODE,I4,
     . 54H NEVER APPEARS IN THE LIST OF NODAL CONNECTION (LNODS))
  925 FORMAT(/,' AND IT HAS NON-ZERO COORDINATES')
      END
