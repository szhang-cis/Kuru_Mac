      SUBROUTINE STRINP
C***********************************************************************
C
C**** THIS ROUTINE READS THE STRATEGY CONTROLLING PARAMETERS FOR THE
C     CURRENT TIME INTERVAL
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'auxl_om.f'
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'inpo_om.f'
C
      PARAMETER (MCOMM=17)
      CHARACTER*5 COMMD(MCOMM)
      DATA COMMD/'ACCEL','ALGOR','ARCLE','CONVE','FINAL','LINES',
     .           'LOADI','OUTPU','DRIFT','NODAL','STEPP','TIMEI',
     .           'RESTA','POSTP','PLOT' ,'RETUR','END_S'/
C
      WRITE(LURES,900)
C
C**** RESET RESTART DEFAULT TO NSTEP
C
      IF(KSAVE.NE.-1) THEN
       NOUTP(1)=NSTEP
       KSAVE=0
      ENDIF
C
C**** READ NEW PARAMETERS IF DESIRED
C
      NPRIN=0
      ITAPE=LUDAT
  200 CALL LISTEN('STRINP',NPRIN,ITAPE)
C
C**** IDENTIFY COMMAND
C
      DO ICOMM=1,MCOMM
       IF(WORDS(1).EQ.COMMD(ICOMM)) GOTO 300
      ENDDO
      GO TO 1000 
C
  300 CONTINUE
      GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), ICOMM
C
C**** 'ACCEL' ACTIVATES OR DISACTIVATES SOLUTION ACCELARATION
C     (DEFAULT LACCE=0)
C
    1 CONTINUE
C
      LACCE=0
      NACCE=INT(PARAM(1))
      IF(WORDS(2).EQ.'ON') LACCE=NACCE
C
      GO TO 200
C
C**** 'ALGOR' DEFINES TYPE OF SOLUTION ALGORITHM 
C     (DEFAULT: NALGO=1)
C
    2 CONTINUE
C
      NALGO=1
      KSTIG=1
      IF(WORDS(2).EQ.'UPDAT') THEN
       IF(INT(PARAM(1)).GT.1) THEN
        NALGO=INT(PARAM(1))
        IF(INT(PARAM(2)).GT.0) KALGO=INT(PARAM(2))
       ENDIF
       IF(WORDS(3).EQ.'NO_GE') KSTIG=0  ! no geometric part of stiffness
      ENDIF
C
      GO TO 200
C
C**** 'ARCLE' SELECTS ARC-LENGTH OPTIONS
C     (DEFAULT KARCL=0, NODIS,NFDIS,NCDIS= NOT DEFINED)
C
    3 CONTINUE
C
      KARCL=0
      IF(WORDS(2).EQ.'NORMA') KARCL=1
      IF(WORDS(2).EQ.'UPDAT') KARCL=2
      IF(WORDS(2).EQ.'SPHER') KARCL=3
      IF(WORDS(2).EQ.'DISPL') KARCL=4
      IF(KARCL.EQ.4) THEN
       NODIS=INT(PARAM(1))
       NFDIS=INT(PARAM(2))
       NCDIS=(NODIS-1)*NDOFN+NFDIS
      ENDIF
C
      GO TO 200
C
C**** 'CONVE' DEFINES CONVERGENCE CRITERIA
C     (DEFAULT KCONV=1, MINTE=50 & TOLER=1%)
C
    4 CONTINUE
C
      KCONV=1
      IF(WORDS(2).EQ.'INCRE') KCONV=2
      IF(WORDS(2).EQ.'DISPL') KCONV=3
      IF(WORDS(2).EQ.'ENERG') KCONV=4
      IF(PARAM(1).GT.0.0D0)   MITER=INT(PARAM(1))
      IF(PARAM(2).NE.0.0D0)   TOLER=PARAM(2)
C
      KCONR=0        ! for problems with zero external & reaction forces
      IF(WORDS(3).EQ.'CONST') KCONR=1
C
      ITRUCM=MITER   ! default: no subincrementation
      TRUCHM=1.0D0
      IF(WORDS(3+KCONR).EQ.'SUBIN') THEN
       IF(NPOIC.GT.0) THEN
        IF(PARAM(3).NE.0.0D0) ITRUCM=INT(PARAM(3))
        IF(PARAM(4).NE.0.0D0) TRUCHM=PARAM(4)
C
        IF(IAUGM.EQ.3) THEN
         IF(PARAM(3).NE.0.0D0) TOLERC=PARAM(3)
        ENDIF
       ENDIF
      ENDIF
C
      GO TO 200
C
C**** 'FINAL' END OF ANALYSIS OPTION
C     (DEFAULT STIFI=0; N.B. IF STIFI=0 IT IS NOT USED AFTERWARDS)
C
    5 CONTINUE
C
      IF(PARAM(1).GT.0.0D0) STIFI=PARAM(1)
      STIFI=STIFI/100.0D0
C
      GO TO 200
C
C**** 'LINES' LINE SEARCH
C
    6 CONTINUE
C
      LINES=0
      IF(WORDS(2).EQ.'ON') LINES=1
C
      GO TO 200
C
C**** 'LOADI' SELECTS LOAD INCREMENTATION OPTIONS
C     (DEFAULT LAUTO=0, DITER= NOT DEFINED)
C
    7 CONTINUE
C
      LAUTO=0
      IF(WORDS(2).EQ.'ITERA') LAUTO=1
      IF(WORDS(2).EQ.'STI1')  LAUTO=2
      IF(WORDS(2).EQ.'STI2')  LAUTO=3
      DITER=PARAM(1)
      IF(LAUTO.GT.1) DITER=DITER/100.0D0
C
      GO TO 200
C
C**** 'OUTPU' SELECTS OUTPUT FREQUENCY
C
    8 CONTINUE
C
      IF(WORDS(2).EQ.'DISPL') THEN   ! + vel. & accel. for dynamic prob.
       NOUTP( 4)=INT(PARAM(1))
       NOUTP(13)=INT(PARAM(2))
      ENDIF
C
      IF(WORDS(2).EQ.'G_STR') THEN
       NOUTP( 5)=INT(PARAM(1))
       NOUTP(14)=INT(PARAM(2))
      ENDIF
      IF(WORDS(2).EQ.'G_PRI') THEN
       NOUTP( 6)=INT(PARAM(1))
       NOUTP(15)=INT(PARAM(2))
      ENDIF
      IF(WORDS(2).EQ.'G_INT') THEN
       NOUTP( 7)=INT(PARAM(1))
       NOUTP(16)=INT(PARAM(2))
      ENDIF
C
      IF(WORDS(2).EQ.'N_STR') THEN
       IF(KSGAU.EQ.0)
     .  CALL RUNMEN('WARNING: NO SMOOTHING CARD HAS BEEN INPUT')
       NOUTP( 8)=INT(PARAM(1))
       NOUTP(17)=INT(PARAM(2))
      ENDIF
      IF(WORDS(2).EQ.'N_PRI') THEN
       IF(KSGAU.EQ.0)
     .  CALL RUNMEN('WARNING: NO SMOOTHING CARD HAS BEEN INPUT')
       NOUTP( 9)=INT(PARAM(1))
       NOUTP(18)=INT(PARAM(2))
      ENDIF
      IF(WORDS(2).EQ.'N_DEF') THEN
       IF(KSGAU.EQ.0)
     .  CALL RUNMEN('WARNING: NO SMOOTHING CARD HAS BEEN INPUT')
       NOUTP(24)=INT(PARAM(1))
       NOUTP(25)=INT(PARAM(2))
      ENDIF
      IF(WORDS(2).EQ.'N_DEP') THEN
       IF(KSGAU.EQ.0)
     .  CALL RUNMEN('WARNING: NO SMOOTHING CARD HAS BEEN INPUT')
       NOUTP(26)=INT(PARAM(1))
       NOUTP(27)=INT(PARAM(2))
      ENDIF
      IF(WORDS(2).EQ.'N_INT') THEN
       IF(KSGAU.EQ.0)
     .  CALL RUNMEN('WARNING: NO SMOOTHING CARD HAS BEEN INPUT')
       NOUTP(20)=INT(PARAM(1))
       NOUTP(21)=INT(PARAM(2))
      ENDIF
C
      IF(WORDS(2).EQ.'REACT') THEN
       NOUTP(22)=INT(PARAM(1))
       NOUTP(23)=INT(PARAM(2))
      ENDIF
C
      IF(WORDS(2).EQ.'N_GAP') THEN
c      IF(ITERME.LE.0) THEN
c       CALL RUNMEN('WARNING: GAP OUTPUT ONLY FOR BIDIRECTIONAL COUPLED
c     .PROBLEMS                     ')
c      ELSE
c       IF(ITERMG.EQ.0)
c     .  CALL RUNMEN('WARNING: GAP OUTPUT ONLY FOR BIDIRECTIONAL COUPLED
c     .PROBLEMS WITH GAP DEPENDENCY ')
c      ENDIF
       NOUTP(28)=INT(PARAM(1))
       NOUTP(29)=INT(PARAM(2))
      ENDIF
C
      IF(WORDS(2).EQ.'HOMOG') THEN
       NOUTP(30)=INT(PARAM(1))
       NOUTP(31)=INT(PARAM(2))
      ENDIF
C
      IF(WORDS(2).EQ.'ALL') THEN
       IF(KSGAU.EQ.0)
     .  CALL RUNMEN('WARNING: NO SMOOTHING CARD HAS BEEN INPUT')
       NOUTP(10)=INT(PARAM(1))
       NOUTP(19)=INT(PARAM(2))
      ENDIF
C
      GO TO 200
C
C**** 'DRIFT' DEFINES TYPE OF PLASTIC BACK PROJECTION 
C     (DEFAULT: NBACK=1)
C     
    9 CONTINUE
C
      NBACK=1
      IF(WORDS(2).EQ.'MODIF') NBACK=2
C
      GO TO 200
C
C**** 'NODAL' ACTIVATES OR DISACTIVATES SMOOTHING OF EXPWP
C     (DEFAULT KSMUS= OFF)
C
   10 CONTINUE
C
      IF(WORDS(2).EQ.'ON')  KSMUS= IABS(KSMUS)
      IF(WORDS(2).EQ.'OFF') KSMUS=-IABS(KSMUS)
      IF(PARAM(1).NE.0.)    WLUMP=PARAM(1)
C
      GO TO 200
C
C**** 'STEPP' DEFINE TIME STEPPING ALGORITHM 
C     (DEFAULT KINTE=0, TALFA & TBETA NOT DEFINED)
C
   11 CONTINUE
C
      KINTE=0
      IF(WORDS(2).EQ.'WILSO') KINTE=1
      IF(WORDS(2).EQ.'NEWMA') KINTE=2
      IF(WORDS(2).EQ.'CENTR') KINTE=3
      IF(WORDS(2).EQ.'HOUBO') KINTE=4         ! not implemented
      IF(WORDS(2).EQ.'HILBE') KINTE=5
      IF(WORDS(2).EQ.'COLLO') KINTE=6
C
      IF(KINTE.EQ.1) THEN
       CONST=PARAM(1)                         ! theta
       IF(CONST.LT.1.37D0) CONST=1.37D0
       TALFA=6.0D0/(CONST*CONST)
       TBETA=3.0D0/CONST
       TGAMA=1.0D0                            ! stiffness in jacobian
      ENDIF
      IF(KINTE.EQ.2) THEN
       CNEW1=PARAM(1)                         ! delta
       IF(CNEW1.LT.0.5D0) CNEW1=0.5D0
       CNEW2=PARAM(2)                         ! alfa
       IF(CNEW2.EQ.0.0D0) CNEW2=(CNEW1+0.5D0)**2/4.0D0
       TALFA=1.0D0/CNEW2                      !     1/alfa
       TBETA=CNEW1/CNEW2                      ! delta/alfa
       TGAMA=1.0D0                            ! stiffness in jacobian
       NCETA=0
       IF(PARAM(3).NE.0.0D0) NCETA=INT(PARAM(3))
      ENDIF
      IF(KINTE.EQ.3) THEN
       TALFA=1.0D0
       TBETA=0.5D0
       TGAMA=0.0D0                            ! no stiffness in jacobian
       NCETA=0
       IF(PARAM(1).NE.0.0D0) NCETA=INT(PARAM(1))
      ENDIF
      IF(KINTE.EQ.4) THEN
       call runend('ERROR: KINTE=4 not implemented')
      ENDIF
      IF(KINTE.EQ.5) THEN
       CNEW3=PARAM(1)                         ! gama (alfa-HHT)
       CNEW1=PARAM(2)                         ! delta
       CNEW2=PARAM(3)                         ! alfa
       IF(CNEW3.GT.0.0D0) THEN
        CNEW1=0.50D0                          ! Newmark
        CNEW2=0.25D0
        CNEW3=0.00D0
       ELSE
        IF(CNEW1.EQ.0.0D0.OR.CNEW2.EQ.0.0D0) THEN
         CNEW1=(1.0D0-2*CNEW3)/2.0D0          ! Hughes' formulas
         CNEW2=(1.0D0-CNEW3)*(1.0D0-CNEW3)/4.0D0
        ENDIF
       ENDIF
       TALFA=1.0D0/CNEW2                      !     1/alfa
       TBETA=CNEW1/CNEW2                      ! delta/alfa
       TGAMA=1.0D0+CNEW3                      ! stiffness in jacobian
       NCETA=0
       IF(PARAM(4).NE.0.0D0) NCETA=INT(PARAM(4))
      ENDIF
      IF(KINTE.EQ.6) THEN
       CNEW3=PARAM(1)                         ! theta
       CNEW1=PARAM(2)                         ! delta
       CNEW2=PARAM(3)                         ! alfa
       IF(CNEW3.LT.1.0D0) THEN
        CNEW1=0.50D0                          ! Newmark
        CNEW2=0.25D0
        CNEW3=1.00D0
       ELSE
        IF(CNEW1.EQ.0.0D0.OR.CNEW2.EQ.0.0D0) THEN
         CNEW1=0.5D0                          ! Hughes' formulas
         CNEW2=CNEW3/(2.0D0*(CNEW3+1.0D0))    ! lower bound
        ENDIF
       ENDIF
       TALFA=1.0D0/(CNEW2*CNEW3*CNEW3)        !     1/(alfa*theta*theta)
       TBETA=CNEW1/(CNEW2*CNEW3)              ! delta/(alfa*theta)
       TGAMA=1.0D0                            ! stiffness in jacobian
       TDELT=CNEW3                            ! theta
       NCETA=0
       IF(PARAM(4).NE.0.0D0) NCETA=INT(PARAM(4))
      ENDIF
      GO TO 200
C
C**** 'TIMEI' SELECTS TIME INCREMENTATION OPTIONS
C     (DEFAULT KOPTI=0, XTIME=1.0D0)
C
   12 CONTINUE
C
      KOPTI=INT(PARAM(1))
      IF(PARAM(2).NE.0.0D0) XTIME=PARAM(2)
      GO TO 200
C
C**** 'RESTA' DEFINES FREQUENCY OF RESTART DUMPING
C     (DEFAULT:  KEEP ONLY LAST STEP )
C
   13 CONTINUE
C
      IF(PARAM(1).NE.0.0D0)  NOUTP( 1)=INT(PARAM(1))
      IF(WORDS(2).EQ.'NONE') KSAVE=-1
      IF(WORDS(2).EQ.'SAVE') KSAVE= 1
      GO TO 200
C
C**** 'POSTP' DEFINES FREQUENCY OF POST-PROCESSING DUMPING
C     (DEFAULT:  NO POST-PROCESSING)
C
   14 CONTINUE
C
C**** ESTABLISHES SOME INDEXES (see setdat.f)
C
      KPPCG=0
      KPPCN=1
C
      IF(KPOST.EQ.1) THEN
       IF(KSGAU.EQ.0)
     .  CALL RUNMEN('WARNING: NO SMOOTHING CARD HAS BEEN INPUT')
       NOUTP( 2)=INT(PARAM(1))
       NOUTP(11)=INT(PARAM(2))
      ENDIF
      GO TO 200
C
C**** 'PLOT' DEFINES X-Y PLOTS
C
   15 NFILE=LUDAT
      IFLAG=0
      IF(WORDS(2).EQ.'NEW') IFLAG=1
      IF(PARAM(1).NE.0.0D0) NFILE=INT(PARAM(1))
      CALL PLOINP(NFILE,IFLAG)
      GO TO 200
C
C**** RETURN MAPPING ALGORITHM
C
   16 IW=2
      IP=1
      IF(WORDS(IW).EQ.'TANGE') THEN   ! tangent form of computing lambda
       IW=IW+1
      ENDIF
      IF(WORDS(IW).EQ.'SECAN') THEN   ! secant form of computing lambda
       IW=IW+1
       NFORZ=2
      ENDIF
      IF(WORDS(IW).EQ.'CONVE') THEN   ! print plastic residual
       IW=IW+1
       IF(PARAM(IP).NE.0.0D0) THEN    ! number of iterations of pl. alg.
        MKONTX=INT(PARAM(IP))
        IF(MKONTX.GT.MKONT)
     .   CALL RUNEND('STRINP: INCREASE MKONT IN setdat.f   ')
        MKONT=MKONTX
        IP=IP+1
       ENDIF
       IF(PARAM(IP).NE.0.0D0) THEN    ! tolerance of pl. alg.
        TOPLA=PARAM(IP)
        IF(TOPLA.LT.0.0D0)
     .   CALL RUNEND('STRINP: WRONG VALUE FOR TOPLA')
        IP=IP+1
       ENDIF
      ENDIF
      IF(WORDS(IW).EQ.'STEPP') THEN   ! stepping strategy for pl. alg.
       IW=IW+1
       ALFAP=PARAM(IP)
       IP=IP+1
      ENDIF
      IF(WORDS(IW).EQ.'SUBIN') THEN   ! subincrementation for pl. alg.
       IW=IW+1
       IF(PARAM(IP).NE.0.0D0) THEN    ! number of subincrements
        MSUBP=INT(PARAM(IP))
        IF(MSUBP.LT.1.0D0)
     .   CALL RUNEND('STRINP: ERROR IN MSUBP VALUE')
        IP=IP+1
       ENDIF
       IF(PARAM(IP).NE.0.0D0) THEN    ! subincrementation factor
        EXCTP=PARAM(IP)
        IF(EXCTP.LT.0.0D0)
     .   CALL RUNEND('STRINP: ERROR IN EXCTP VALUE')
        IP=IP+1
       ENDIF
      ENDIF
      IF(WORDS(IW).EQ.'ZERO_') THEN   ! zero iteration computation
       IW=IW+1
       IF(PARAM(IP).NE.0.0D0) THEN    ! indicator
        ITEPL=INT(PARAM(IP))
        IF(ITEPL.LT.-1.0D0.OR.ITEPL.GT.0.0D0)
     .   CALL RUNEND('STRINP: ERROR IN ITEPL VALUE')
        IP=IP+1
       ENDIF
      ENDIF
      IF(WORDS(IW).EQ.'CONSI') THEN   ! consistent constitutive tensor
       IW=IW+1
       IF(PARAM(IP).NE.0.0D0) THEN    ! indicator
        ICOCO=INT(PARAM(IP))
        IF(ICOCO.LT.-1.0D0.OR.ICOCO.GT.1.0D0)
     .   CALL RUNEND('STRINP: ERROR IN ICOCO VALUE')
        IP=IP+1
       ENDIF
      ENDIF
      IF(WORDS(IW).EQ.'NO_UP') THEN   ! does not update Cep
       IW=IW+1
       NALGP=0
      ENDIF
      GO TO 200
C
C**** END
C
   17 CONTINUE
C
      RETURN
 1000 CALL RUNEND('STRINP:ERROR IN STRATEGY DATA BLOCK')
  900 FORMAT(1H1,///)
      END
