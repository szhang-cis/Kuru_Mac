      SUBROUTINE INPCEK(COORD,LNODS,MATNO,NDIME,NELEM,NMATS,NNODE,
     .                  NPOIN,PROEL,NPREL,NGRUP,NPOIC,NNODC,NOCOI,LURES,
     .                  KGAPC,LARGE,LARGC)
C***********************************************************************
C
C**** THIS ROUTINE CHECKS THE GEOMETRICAL INPUT DATA
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION COORD(NDIME,*), LNODS(NNODE,*),
     .          MATNO(*)
      DIMENSION NEROR(24),      PROEL(NPREL,NGRUP)
C
C**** COMPUTES NUMBER OF POINTS WITHOUT CONTACT (AL METHOD)
C
      NPOI1=NPOIN-NPOIC
C
C**** INITIALISED THE LIST OF ERRORS
C
      DO 5 IEROR=13,24
    5 NEROR(IEROR)=0
C
C**** REDEFINES LOCAL NODE NUMBER OF CONTACT ELEMENTS FOR NON-COINCIDENT
C     MESH & REASSIGNS CONTACT PARAMETERS
C
C     Notes:
C
C     PROEL(19,IGRUP): contact set (ISETC)
C     PROEL(20,IGRUP): contact index (IMAES)
C                           (1=no self-contact, 2=self-contact & 3=both)
C     These two indexes are useful to the contact searching operations
C
C     PROEL(21,IGRUP): real NNODL (NNODN=NNODL+NGAUL*NNODL)
C
      IF(NOCOI.EQ.1) THEN             ! non coincident mesh
       DO IGRUP=1,NGRUP
        NNODL=INT(PROEL(2,IGRUP))
        ITYPE=INT(PROEL(5,IGRUP))
        IF(ITYPE.EQ.32) THEN
         IF(INT(PROEL(19,IGRUP)).EQ.0.AND.
     .      INT(PROEL(20,IGRUP)).EQ.0) THEN             ! old version
          PROEL(2,IGRUP)=FLOAT(NNODL-2)                 ! new NNODL
         ENDIF
        ENDIF
       ENDDO
       DO IELEM=1,NELEM
        IGRUP=MATNO(IELEM)
        ITYPE=INT(PROEL(5,IGRUP))
        IF(ITYPE.EQ.32) THEN
         NNODL=INT(PROEL(2,IGRUP))
         IF(INT(PROEL(19,IGRUP)).EQ.0.AND.
     .      INT(PROEL(20,IGRUP)).EQ.0) THEN             ! old version
          PROEL(19,IGRUP)=FLOAT(LNODS(NNODL+1,IELEM))   ! contact set
          PROEL(20,IGRUP)=FLOAT(LNODS(NNODL+2,IELEM))   ! contact index
         ENDIF
C
         IF(NDIME.EQ.2) THEN                            ! some checks
          IF(NNODL.GT.2)
     .     CALL RUNEND('ERROR: NNODL GT 2 FOR ITYPE=32')
         ENDIF
         IF(NDIME.EQ.3) THEN
          IF(NNODL.GT.4)
     .     CALL RUNEND('ERROR: NNODL GT 4 FOR ITYPE=32')
         ENDIF
C
         NGAUL=INT(PROEL(4,IGRUP))
         PROEL(21,IGRUP)=FLOAT(NNODL+NGAUL*NNODL)       ! real NNODL
         NNODN=INT(PROEL(21,IGRUP))
         IF(NNODN.GT.NNODE)
     .    CALL RUNEND('ERROR: NNODN GT NNODE => CHANGE INPUT DATA')
C
         DO INODL=NNODL+1,NNODN
          LNODS(INODL,IELEM)=1        ! proper initialization for solver
         ENDDO
        ENDIF
       ENDDO
      ENDIF
C
      IF(NOCOI.EQ.2) THEN             ! coincident & non coincident mesh
       DO IELEM=1,NELEM
        IGRUP=MATNO(IELEM)
        ITYPE=INT(PROEL(5,IGRUP))
        IF(ITYPE.EQ.32) THEN
         NOCOL=INT(PROEL(22,IGRUP))
         IF(NOCOL.EQ.1) THEN          ! non coincident mesh
          NNODL=INT(PROEL(2,IGRUP))
C
          IF(NDIME.EQ.2) THEN                           ! some checks
           IF(NNODL.GT.2)
     .      CALL RUNEND('ERROR: NNODL GT 2 FOR ITYPE=32')
          ENDIF
          IF(NDIME.EQ.3) THEN
           IF(NNODL.GT.4)
     .      CALL RUNEND('ERROR: NNODL GT 4 FOR ITYPE=32')
          ENDIF
C
          NGAUL=INT(PROEL(4,IGRUP))
          PROEL(21,IGRUP)=FLOAT(NNODL+NGAUL*NNODL)      ! real NNODL
          NNODN=INT(PROEL(21,IGRUP))
          IF(NNODN.GT.NNODE)
     .     CALL RUNEND('ERROR: NNODN GT NNODE => CHANGE INPUT DATA')
C
          DO INODL=NNODL+1,NNODN
           LNODS(INODL,IELEM)=1       ! proper initialization for solver
          ENDDO
         ENDIF
        ENDIF
       ENDDO
      ENDIF
C
C**** CHECK ORDER OF ELEMENT DEFINITION WHEN A DISCONTINUOUS GALERKIN
C     (DG) FORMULATION IS USED (DG ELEMENTS MUST BE DEFINED AFTER SOLID
C     ELEMENTS
C
      IF(IGALE.EQ.1) THEN
       IG=0
       DO IELEM=1,NELEM
        IGRUP=MATNO(IELEM)
        ITYPE=INT(PROEL(5,IGRUP))
        IF(ITYPE.EQ.30.AND.IG.EQ.1)
     .   CALL RUNEND('ERROR IN ORDER OF ELEMENT DEFINITION FOR DG')
        IF(ITYPE.EQ.33) IG=1
       ENDDO
      ENDIF
C
C**** CHECK AGAINST TWO IDENTICAL NONZERO NODAL COORDINATES
C
      TOLERG=0.0D0          ! better as input
      DO 40 IPOIN=2,NPOI1
C$DIR SCALAR
      DO 30 JPOIN=1,IPOIN-1
C$DIR SCALAR
      DO 20 IDIME=1,NDIME
c  20 IF(COORD(IDIME,IPOIN).NE.COORD(IDIME,JPOIN)) GOTO 30
   20 IF(ABS(COORD(IDIME,IPOIN)-COORD(IDIME,JPOIN)).GT.TOLERG) GOTO 30
C
C**** FORM TO AVOID CHECK OF COINCIDENT COORDINATES FOR MECHANICAL
C     CONTACT ELEMENTS
C
      IF(NOCOI.EQ.0.OR.NOCOI.EQ.2) THEN   ! coincident mesh
       jjj=0
       do ielem=1,nelem
        igrup=matno(ielem)
        itype=int(proel(5,igrup))
        NNODL=INT(PROEL(2,IGRUP))
C
        if(itype.eq.3.or.itype.eq.4) then
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
          CALL RUNEND('INPCEK:ERROR DETECTED WITH NNODL.NE.2 ')
         endif
        endif
C
        if(itype.eq.32.or.itype.eq.33) then
         NOCOL=INT(PROEL(22,IGRUP))
         IF(NOCOL.EQ.0) THEN
          if(ndime.eq.1) then
           if(nnodl.eq.2) then
            iii=0
            do inode=1,nnode
             if(lnods(inode,ielem).eq.ipoin) iii=iii+1
             if(lnods(inode,ielem).eq.jpoin) iii=iii+1
            enddo
            if(iii.eq.2) go to 30   ! same element
            if(iii.eq.1) then       ! different element
             jjj=jjj+1
             if(jjj.eq.2) go to 30
            endif
           else
            CALL RUNEND('INPCEK:ERROR DETECTED WITH NNODL.NE.2 ')
           endif
          endif             ! ndime.eq.1
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
            if(iii.eq.2) then       ! same element
             if(inodi.eq.0.or.inodj.eq.0)
     .        call runend('inodi-inodj = 0       ')
             if(inodi.eq.1.and.inodj.eq.4) go to 30
             if(inodi.eq.4.and.inodj.eq.1) go to 30
             if(inodi.eq.2.and.inodj.eq.3) go to 30
             if(inodi.eq.3.and.inodj.eq.2) go to 30
            endif
            if(iii.eq.1) then       ! different element
             jjj=jjj+1
             if(jjj.eq.2) go to 30
            endif
           else
            CALL RUNEND('INPCEK:ERROR DETECTED WITH NNODL.NE.4 ')
           endif
          endif             ! ndime.eq.2
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
            if(iii.eq.2) then       ! same element
             if(inodi.eq.0.or.inodj.eq.0)
     .        call runend('inodi-inodj = 0       ')
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
            if(iii.eq.1) then       ! different element
             jjj=jjj+1
             if(jjj.eq.2) go to 30
            endif
           else
            CALL RUNEND('INPCEK:ERROR DETECTED WITH NNODL.NE.6-8 ')
           endif
          endif             ! ndime.eq.3
         ENDIF              ! NOCOL=0
        endif               ! itype.eq.32.or.itype.eq.33
       enddo                ! ielem=1,nelem
      ENDIF                 ! NOCOI=0 OR NOCOI=2
C
      IF(NOCOI.EQ.1.OR.NOCOI.EQ.2) THEN   ! non coincident mesh
       ielem1=0
       igrup1=0
       isetc1=0
       imaes1=0
       do ielem=1,nelem
        igrup=matno(ielem)
        itype=int(proel(5,igrup))
        if(itype.eq.32) then
         NOCOL=INT(PROEL(22,IGRUP))
         IF(NOCOL.EQ.1) THEN
          NNODL=INT(PROEL(2,IGRUP))
          isetc=int(proel(19,igrup))
          imaes=int(proel(20,igrup))
          do inodl=1,nnodl
           if(lnods(inodl,ielem).eq.ipoin) then
            ielem1=ielem
            igrup1=igrup
            isetc1=isetc
            imaes1=imaes
           endif
          enddo
         ENDIF
        endif
       enddo
C
       if(ielem1.ne.0) then
        ielem2=0
        igrup2=0
        isetc2=0
        imaes2=0
        do ielem=1,nelem
         if(ielem.ne.ielem1) then
          igrup=matno(ielem)
          itype=int(proel(5,igrup))
          if(itype.eq.32) then
           NOCOL=INT(PROEL(22,IGRUP))
           IF(NOCOL.EQ.1) THEN
            NNODL=INT(PROEL(2,IGRUP))
            isetc=int(proel(19,igrup))
            imaes=int(proel(20,igrup))
            do inodl=1,nnodl
             if(lnods(inodl,ielem).eq.jpoin) then
              ielem2=ielem
              igrup2=igrup
              isetc2=isetc
              imaes2=imaes
             endif
            enddo
           ENDIF
          endif
         endif
        enddo
C
        if(ielem2.ne.0) then
         if(isetc1.eq.isetc2) then
          if(igrup1.eq.igrup2) then
           if((imaes1.eq.2.and.imaes2.eq.2).or.
     .        (imaes1.eq.3.and.imaes2.eq.3)) go to 30
          else
           if(imaes1.eq.1.and.imaes2.eq.1) go to 30
          endif
         endif
        endif
       endif
      ENDIF                 ! NOCOI=1 OR NOCOI=2
C
      NEROR(13)=NEROR(13)+1
      WRITE(LURES,900) IPOIN,JPOIN
   30 CONTINUE
   40 CONTINUE
C
C**** COMPATIBLE INPUT WITH THE REST OF ELEMENTS (ONLY FOR 2D COINCIDENT
C     MESH AT CONTACT)
C
C     INPUT 4-----3         INTERNAL 3-----4
C           .     .                  .     .
C           1-----2                  1-----2
C
      IF(NOCOI.EQ.0.OR.NOCOI.EQ.2) THEN   ! coincident mesh
       if(ndime.eq.2) then
        do ielem=1,nelem
         igrup=matno(ielem)
         itype=int(proel(5,igrup))
         NNODL=INT(PROEL(2,IGRUP))
         if(itype.eq.32.or.itype.eq.33) then
          NOCOL=INT(PROEL(22,IGRUP))
          IF(NOCOL.EQ.0) THEN
           if(nnodl.eq.4) then
            lauxi=lnods(3,ielem)
            lnods(3,ielem)=lnods(4,ielem)
            lnods(4,ielem)=lauxi
           else
            CALL RUNEND('INPCEK:ERROR DETECTED WITH NNODL.NE.4 ')
           endif
          ENDIF
         endif
        enddo                ! ielem=1,nelem
       endif                 ! ndime.eq.2
      ENDIF
C
C**** MODIFIES "LNODS" FOR AUMENTED LAGRANGIAN APPROACH
C
      if(npoic.ne.0) then
       iii=npoi1
       do ielem=1,nelem
        igrup=matno(ielem)
        itype=int(proel(5,igrup))
        NNODL=INT(PROEL(2,IGRUP))
C
        if(itype.eq.4) then
         iii=iii+1
         PROEL(2,IGRUP)=3.0          ! nnodl
         lnods(3,ielem)=iii
         if(nnode.lt.3)
     .    CALL RUNEND('INPCEK: ERROR NNODE CAN NOT BE LT 3 FOR ALM')
        endif
C
        if(itype.eq.32) then
         if(ndime.eq.1) then
          IF(NNODC.NE.1)
     .     CALL RUNEND('INPCEK: ERROR IN NNODC')
          iii=iii+1
          lnods(3,ielem)=iii
          if(nnode.lt.3)
     .     CALL RUNEND('INPCEK: ERROR NNODE CAN NOT BE LT 3 FOR ALM')
         endif
         if(ndime.eq.2) then
          if(nnodl.eq.4) then
           IF(NNODC.NE.2)
     .      CALL RUNEND('INPCEK: ERROR IN NNODC')
C
C**** ASSIGNS TWO NEW POINTS PER CONTACT ELEMENT
C
           iii=iii+1
           lnods(5,ielem)=iii
           iii=iii+1
           lnods(6,ielem)=iii
C
           iccc=1
           if(iccc.eq.0) then              ! no funca

           inod1=lnods(1,ielem)
           inod2=lnods(2,ielem)
           inod3=lnods(3,ielem)
           inod4=lnods(4,ielem)
           i56=4
           do jelem=1,ielem-1
            jgrup=matno(jelem)
            jtype=int(proel(5,jgrup))
            if(jtype.eq.32) then
             jnod1=lnods(1,jelem)
             jnod2=lnods(2,jelem)
             jnod3=lnods(3,jelem)
             jnod4=lnods(4,jelem)
             if(inod1.eq.jnod1.and.inod3.eq.jnod3) then
              i56=i56+1
              lnods(5,ielem)=lnods(5,jelem)
             endif
             if(inod1.eq.jnod3.and.inod3.eq.jnod1) then
              i56=i56+1
              lnods(5,ielem)=lnods(5,jelem)
             endif
             if(inod1.eq.jnod2.and.inod3.eq.jnod4) then
              i56=i56+1
              lnods(5,ielem)=lnods(6,jelem)
             endif
             if(inod1.eq.jnod4.and.inod3.eq.jnod2) then
              i56=i56+1
              lnods(5,ielem)=lnods(6,jelem)
             endif
             if(inod2.eq.jnod2.and.inod4.eq.jnod4) then
              i56=i56+2
              lnods(6,ielem)=lnods(6,jelem)
             endif
             if(inod2.eq.jnod4.and.inod4.eq.jnod2) then
              i56=i56+2
              lnods(6,ielem)=lnods(6,jelem)
             endif
             if(inod2.eq.jnod1.and.inod4.eq.jnod3) then
              i56=i56+2
              lnods(6,ielem)=lnods(5,jelem)
             endif
             if(inod2.eq.jnod3.and.inod4.eq.jnod1) then
              i56=i56+2
              lnods(6,ielem)=lnods(5,jelem)
             endif
            endif
           enddo
           if(i56.eq.4) then
            iii=iii+1
            lnods(5,ielem)=iii
            iii=iii+1
            lnods(6,ielem)=iii
           endif
           if(i56.eq.5) then
            iii=iii+1
            lnods(6,ielem)=iii
           endif
           if(i56.eq.6) then
            iii=iii+1
            lnods(5,ielem)=iii
           endif

           endif            ! iccc
c
           if(nnode.lt.6)
     .      CALL RUNEND('INPCEK: ERROR NNODE CAN NOT BE LT 3 FOR ALM')
          endif
          if(nnodl.eq.6) then
           CALL RUNEND('INPCEK: ERROR NNODE=6')
          endif
         endif              ! ndime.eq.2
         if(ndime.eq.3) then
          if(nnodl.eq.6) then
           CALL RUNEND('INPCEK: ERROR NNODE=6')
          endif
          if(nnodl.eq.8) then
           CALL RUNEND('INPCEK: ERROR NNODE=8')
          endif
          if(nnodl.eq.16) then
           CALL RUNEND('INPCEK: ERROR NNODE=16')
          endif
          if(nnodl.eq.18) then
           CALL RUNEND('INPCEK: ERROR NNODE=18')
          endif
         endif              ! ndime.eq.3
        endif               ! itype.eq.32
       enddo                ! ielem=1,nelem
C
C**** RECOMPUTES "NNODL"
C
       do igrup=1,ngrup
        itype=int(proel(5,igrup))
        if(itype.eq.32) then
         nnodl=int(proel(2,igrup))
         if(ndime.eq.1) then
          if(nnodl.eq.2) proel(2,igrup)=3.0        ! nnodl
         endif
         if(ndime.eq.2) then
          if(nnodl.eq.4) proel(2,igrup)=6.0        ! nnodl
          if(nnodl.eq.6) proel(2,igrup)=9.0        ! nnodl
         endif
         if(ndime.eq.3) then
          if(nnodl.eq.6) proel(2,igrup)=9.0        ! nnodl
          if(nnodl.eq.8) proel(2,igrup)=12.0       ! nnodl
         endif
        endif
       enddo
C
       if(iii.gt.npoi1) then
        jjj=npoi1+npoic
        if(iii.ne.jjj)
     .   CALL RUNEND('INPCEK:ERROR DETECTED WITH NPOIC')
       else
        npoin=npoin-npoic   ! npoic ne 0 but no contact elements
       endif
      endif                 ! npoic.ne.0
C
C**** CONTROLS THE EXISTENCE OF GAP ELEMENTS FOR OUTPUT & POSTPROCESS
C
      KGAPC=0
      IF(NPOIC.NE.0.OR.NOCOI.GT.0) THEN
       KGAPC=1
      ELSE
       DO IELEM=1,NELEM
        LGRUP=MATNO(IELEM)
        LTYPE=INT(PROEL(5,LGRUP))
        IF(LTYPE.EQ.4.OR.LTYPE.EQ.32) KGAPC=1
       ENDDO
      ENDIF
C
      IF(KGAPC.EQ.1) THEN
       IF(NOCOI.EQ.0.OR.NOCOI.EQ.2) THEN
        IF(LARGE.NE.0)
     .   CALL RUNMEN('WARNING: COINCIDENT CONTACT MESH WITH LARGE NE 0')
        IF(LARGC.NE.0)
     .   CALL RUNMEN('WARNING: COINCIDENT CONTACT MESH WITH LARGC NE 0')
       ENDIF
       IF(NOCOI.EQ.1.OR.NOCOI.EQ.2) THEN
        LARGX=0
        IF(LARGE.NE.0) LARGX=1
        IF(LARGX.NE.LARGC)
     .   CALL RUNMEN('WARNING: LARGE NE LARGC')
       ENDIF
      ENDIF
C
C**** CHECK THE LIST OF ELEMENT PROPERTY NUMBERS
C
C$DIR SCALAR
      DO 50 IELEM=1,NELEM
      IF((MATNO(IELEM).LE.0).OR.(MATNO(IELEM).GT.NGRUP)) THEN
       NEROR(14)=NEROR(14)+1
       WRITE(LURES,905) IELEM,MATNO(IELEM)
      ENDIF
   50 CONTINUE
C
C**** CHECK FOR IMPOSSIBLE NODE NUMBERS
C
      DO 60 IELEM=1,NELEM
C$DIR SCALAR
      IGRUP=MATNO(IELEM)
      NNODL=INT(PROEL(2,IGRUP))
      DO 60 INODE=1,NNODL
      IF((LNODS(INODE,IELEM).LT.0).OR.
     .   (LNODS(INODE,IELEM).GT.NPOIN)) THEN
       NEROR(16)=NEROR(16)+1
       WRITE(LURES,910) IELEM,INODE,LNODS(INODE,IELEM),NPOIN
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
C     CALCULATE INCREASE OR DECREASE IN FRONTWIDTH AT EACH ELEMENT STAGE
C
      KSTAR=0
      DO 120 IELEM=1,NELEM
      KZERO=0
      IGRUP=MATNO(IELEM)
      NNODL=INT(PROEL(2,IGRUP))
      DO 130 INODE=1,NNODL
      IF(LNODS(INODE,IELEM).NE.IPOIN) GOTO 130
      KZERO=KZERO+1
C
      IF(KZERO.GT.1) THEN
       NEROR(17)=NEROR(17)+1
       WRITE(LURES,915) IPOIN,IELEM
      ENDIF
C
      IF(KSTAR.EQ.0) KSTAR=IELEM
C
      KLAST=IELEM
      NLAST=INODE
  130 CONTINUE
  120 CONTINUE
C
C**** CHECK IF THIS IS AN UNUSED NODE & IF IT HAS NON-ZERO COORDINATES
C
      IF(KSTAR.EQ.0)THEN
       WRITE(LURES,920) IPOIN
       NEROR(18)=NEROR(18)+1
C
       SIGMA=0.0D0
C$DIR SCALAR
       DO 140 IDIME=1,NDIME
  140  SIGMA=SIGMA+COORD(IDIME,IPOIN)*COORD(IDIME,IPOIN)
       IF(SIGMA.NE.0.)THEN
        NEROR(19)=NEROR(19)+1
        WRITE(LURES,925)
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
     . CALL RUNEND('ERROR: WRONG R COORDINATES')
C
C**** VERIFIED IF ANY INPUT ERRORS HAVE BEEN DETECTED
C
      KEROR=0
      DO 200 IEROR=13,24
  200 IF(NEROR(IEROR).NE.0) KEROR=1
      IF(KEROR.EQ.1)  
     . CALL RUNEND('ERROR DETECTED IN INPCEK           ')
C
      RETURN
  900 FORMAT(1X,
     .'IDENTICAL COORDINATES HAS BEEN FOUND FOR NODAL POINTS N0.',2I5)
  905 FORMAT(1X,'SET NO. OUT OF RANGE 1-NGRUP. MATNO(',I5,')=',I5) 
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
