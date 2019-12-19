      SUBROUTINE PARASSM(CSTIF,VALUES,LCOLUMNS,IRWNDX,LNUEQ,KSYMM,
     .                  LNODS,NDOFN,NELEM,NEQNS,NEVAB,NLAST,NNODE,
     .                  NWIDT,NTOTV,NPOIN,
     .                  NPREL,NGRUP,NPROP,NMATS,
     .                  NDATA,NPREV,NSTAT,NMATX,NDIME,
     .                  NDISR,NDISO,
     .                  MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .                  ELVAR,ELMAT,WORK1,TEMPN,DISPR,DTEMP,
     .                  VNORM,DISTO,COORD,INFRI,COFRI,LACTI,
     .                  NPRDSSIZE,IPRDSINIT,CHANGEPHASE)
C***********************************************************************
C
C**** THIS ROUTINE ASSEMBLES THE ELEMENTAL STIFFNESS MATRIX 
C     (SYMMETRIC OR UNSYMMETRIC) INTO THE GLOBAL STIFFNESS ARRAYS
C
C.... INPUT PARAMETERS
C
C     IRWNDX(NEQNS+1)         - ARRAY WITH INDEX TO FIRST ELEMENT OF THE ROW
C     LCOLUMNS(UNKNOW)        - ARRAY WITH INDEX TO ELEMENTS IN THE ROW  
C     VALUES(UNKNOW)          - VALUES OF THE MATRIX
C     NPRDSSIZE               - APPROX OF NUMBER NON-ZEROS OF THE MATRIX
C
C     KSYMM                   - FLAG DEFINING SYMMETRY OF MATRIX
C     NWIDT                   - MAXIMUM HALF-WIDTH ALLOWED (NO USADO)
C
C**** OUPUT PARAMETERS
C
C     IRWNDX(NEQNS+1)            - DIAGONAL OF GLOBAL MATRIX
C     LCOLUMNS(NPRDSSIZE)     - LOWER TRIANGULAR PART OF PROFILE STORED
C                               GLOBAL MATRIX (NEEDED IF UNSYMMETRIC)
C     VALUES(NPRDSSIZE)     - UPPER TRIANGULAR PART OF PROFILE STORED
C                               GLOBAL MATRIX
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8, ALLOCATABLE, INTENT(INOUT) :: VALUES(:)
      REAL*8, ALLOCATABLE :: TEMP1(:)
      INTEGER*4, ALLOCATABLE, INTENT(INOUT) :: IRWNDX(:)
      INTEGER*4, ALLOCATABLE, INTENT(INOUT) :: LCOLUMNS(:)
      INTEGER*4, ALLOCATABLE :: TEMP2(:)
      INTEGER*4 :: CHANGEPHASE
      
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'auxl_om.f'
      INCLUDE 'inte_om.f'
C
      DIMENSION CSTIF(NEVAB,NEVAB), LNODS(NNODE,*), LNUEQ(*)
C
      DIMENSION MATNO(NELEM),
     .          PROEL(NPREL,NGRUP), PROPS(NPROP,NMATS),
     .          ELDAT(NDATA),       ELPRE(NPREV),
     .          ELVAR(NSTAT),       ELMAT(NMATX),
     .          WORK1(*)
      DIMENSION TEMPN(NPOIN,2),     DISPR(NTOTV,NDISR),
     .          DTEMP(NPOIN),       VNORM(NTOTV)
      DIMENSION DISTO(NTOTV,NDISO), COORD(NDIME,NPOIN)
      DIMENSION INFRI(*),           COFRI(*),
     .          LACTI(NELEM)
     

      INTEGER*4      NOCOI,NNODN,LARGC,LICOI,NSKIC,NOCOL
      COMMON/CONTNCM/NOCOI,NNODN,LARGC,LICOI,NSKIC,NOCOL

      COMMON/LOGUN/LUDTS,LUSOL,LUFRO,LUFRH,LUDAT,LUPRI,LURES,
     .             LUSO2,LUFR2,LUPOS,LURST,LUBFG,LUPIP,LUPAN,LUINF,
     .             LUGEO,LUSET,LUMAT,LUINI,LULOA,LUFIX,LUIN1,
     .             LUTUN,LUCON,LUACT,LUFAN,
     .             LUCU1,LUCU2,LUCU3,LUCU4,LUCU5,LUCU6,LUCU7,
     .             LUCU8,LUCU9,LUC10
    
     
C      if ((IPRDSINIT.eq.0).or.(NOCOI.EQ.1)) then
      CHANGEPHASE=0
      IF ((IPRDSINIT.EQ.0).OR.(NEWBO.EQ.1).OR.(IREST.NE.0)
     .    .OR.(NOCOI.GT.0.AND.LICOI.EQ.1)
     .    .OR.(NOCOI.GT.0.AND.LICOI.EQ.0.AND.MOD(ISTEP,NSKIC).EQ.0))
     . THEN
       CHANGEPHASE=1 
       WRITE(LURES,*) "Rearmando Matriz Dispersa" 
C
C**** INITIALISE ROW & COLUMNS POINTERS
C
 1200     DO IEQNS=1,NEQNS
       IRWNDX(IEQNS)=0
      END DO
      NCOLMAX=NPRDSSIZE/NEQNS
C
C***COMPUTING IRWNDX & LCOLUMNS
C
      DO 1000 IELEM=1,NELEM
      DO 100 INODE=1,NNODE
      IPOIN=LNODS(INODE,IELEM)
      IF(IPOIN.EQ.0) GO TO 1000
      IEVAB=(INODE-1)*NDOFN
      ITOTV=(IPOIN-1)*NDOFN
C
      DO 100 IDOFN=1,NDOFN
      IEVAB=IEVAB+1
      ITOTV=ITOTV+1
      IDEST=LNUEQ(ITOTV)
      IF(IDEST.GT.0) THEN
C
C**** LOOP THROUGH THE COLUMNS TO PERFORM THE ASSEMBLY
C
       DO 40 JNODE=1,NNODE
       JPOIN=LNODS(JNODE,IELEM)
       IF(JPOIN.EQ.0) GO TO 200
       JEVAB=(JNODE-1)*NDOFN
       JTOTV=(JPOIN-1)*NDOFN
       DO 40 JDOFN=1,NDOFN
       JEVAB=JEVAB+1
       JTOTV=JTOTV+1
C----------------------------------------------------------------
       JDEST=LNUEQ(JTOTV)
       IF(JDEST.GT.0) THEN
         K1= (IDEST-1)*NCOLMAX+1
         K2= K1+IRWNDX(IDEST)-1
         DO K=K1,K2
           IF(LCOLUMNS(K).ge.JDEST) GO TO 43
         END DO
   43    IF (K.GT.K2.OR.LCOLUMNS(K).NE.JDEST)  THEN
         IF (KSYMM.EQ.0.OR.JDEST.GE.IDEST) THEN
          IRWNDX(IDEST)=IRWNDX(IDEST)+1
          IF(IRWNDX(IDEST).GT.NCOLMAX) THEN
C***************************************************************
C                                                           Redimensionando values
           LARGE=SIZE(VALUES)
           NPRDSSIZE=INT(LARGE*1.5)
           WRITE(LURES,*) "REALLOCATE VALUES (from, to)",
     .                 LARGE,NPRDSSIZE,';',
     .                 LARGE/(1024.0*1024.0)*8.0,'MB',
     .                 NPRDSSIZE/(1024.0*1024.0)*8.0,'MB'
           DEALLOCATE(VALUES)
           ALLOCATE(VALUES(NPRDSSIZE),STAT=IERROR)
           IF (IERROR.NE.0) 
     .      CALL RUNEND("I CAN'T REALLOCATE VALUES IN PARASSM")
C                                                           Redimensionando lcolumn
          
           WRITE(LURES,*) "REALLOCATE LCOLUMNS (from, to)",
     .                 LARGE,NPRDSSIZE,';',
     .                 LARGE/(1024.0*1024.0)*8.0,'MB',
     .                 NPRDSSIZE/(1024.0*1024.0)*8.0,'MB'
           DEALLOCATE(LCOLUMNS)
           ALLOCATE(LCOLUMNS(NPRDSSIZE),STAT=IERROR)
           IF (IERROR.NE.0) 
     .      CALL RUNEND("I CAN'T REALLOCATE VALUES IN PARASSM")
           
           GOTO 1200
          END IF
          DO KK=K2,K,-1
          LCOLUMNS(KK+1)=LCOLUMNS(KK)
          END DO
          LCOLUMNS(K)=JDEST
         END IF
         END IF
       END IF     
   40  CONTINUE
  200  CONTINUE
      ENDIF
C
  100 CONTINUE
C
 1000 CONTINUE
C
C-------------------------------------------------------------
      do kcoef=1,neqns
          K1= IRWNDX(KCOEF+1)
          IRWNDX(KCOEF+1)=IRWNDX(KCOEF)+K1
      end do
      do kcoef=neqns,1,-1
          IRWNDX(KCOEF+1)=IRWNDX(KCOEF)+1
      end do
      IRWNDX(1)=1


      do kcoef=1,neqns
          K1OLD= (KCOEF-1)*NCOLMAX+1
          KEQNS= IRWNDX(KCOEF+1)-IRWNDX(KCOEF)
          K1NEW= IRWNDX(KCOEF)
          do KMOVE=1,KEQNS
            LCOLUMNS(K1NEW+KMOVE-1)=LCOLUMNS(K1OLD+KMOVE-1)
          end do
      end do
      end if                         !iprdsinit=0

      ncoef=IRWNDX(NEQNS+1)-1
c     write(6,*)'ncoef nuevo=',ncoef
      do kcoef=1,ncoef
        values(kcoef)=0.d0
      end do
C----------------------------------------------



C
C**** ASSEMBLE GLOBAL MATRIX 
C
      DO 2000 IELEM=1,NELEM
C
C**** READ CSTIF FROM DATA BASE
C 
      IF(NMEMO7M.EQ.0) THEN
       CALL DATBAS(CSTIF,    6,    2)
      ELSE
       CALL STIFMS(ELDAT,ELPRE,ELVAR,ELMAT,LNODS,MATNO,PROEL,
     .             PROPS,WORK1,TEMPN(1,1),DISPR(1,1),DTEMP,VNORM,COORD,
     .             DISTO(1,1),INFRI,COFRI,LACTI,CSTIF)
      ENDIF
C
      DO 220 INODE=1,NNODE
      IPOIN=LNODS(INODE,IELEM)
      IF(IPOIN.EQ.0) GO TO 2000
      IEVAB=(INODE-1)*NDOFN
      ITOTV=(IPOIN-1)*NDOFN
C
      DO 220 IDOFN=1,NDOFN
      IEVAB=IEVAB+1
      ITOTV=ITOTV+1
      IDEST=LNUEQ(ITOTV)
      IF(IDEST.GT.0) THEN
        KINIT=IRWNDX(IDEST)
        KEND0=IRWNDX(IDEST+1)-1
C
C**** LOOP THROUGH THE COLUMNS TO PERFORM THE ASSEMBLY
C
       DO 80 JNODE=1,NNODE
       JPOIN=LNODS(JNODE,IELEM)
       IF(JPOIN.EQ.0) GO TO 400
       JEVAB=(JNODE-1)*NDOFN
       JTOTV=(JPOIN-1)*NDOFN
       DO 80 JDOFN=1,NDOFN
       JEVAB=JEVAB+1
       JTOTV=JTOTV+1
C-------------------------------------------------------
       JDEST=LNUEQ(JTOTV)
       IF(JDEST.GT.0) THEN
       IF (KSYMM.EQ.0.OR.JDEST.GE.IDEST) THEN
          K1 = KINIT
          K2 = KEND0+1
 82       if(k2-k1.eq.1) goto 83
          KMIDDLE = (k1+k2)/2
          if (jdest.ge.lcolumns(kmiddle)) then
             k1 = kmiddle
          else
             k2 = kmiddle
          end if
          goto 82
 83    if (lcolumns(k1).ne.jdest) then
          call runend('couldnt find jdest')
       endif
          values(k1)=values(k1)+cstif(ievab,jevab)
       END IF
       ENDIF
   80  CONTINUE
  400  CONTINUE
      ENDIF
C
  220 CONTINUE
C
 2000 CONTINUE
C
C----------------------------------------------------------
C      CALL CHECKSM (VALUES,NCOEF,COEFSUM)
C      write (6,*) 'Para kprob= ',kprob
C      write (6,*) 'ncoef, loc-values in skyasspm, coefsum: ',ncoef,
C     .      loc(values),coefsum
C---------------------------------------------------------


      RETURN
      END
