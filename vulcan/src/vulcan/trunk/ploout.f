      SUBROUTINE PLOOUT(DISTO,ELVAR,TLOAD,eldat,matno,proel,props,REFOR,
     .                  CARTD,ELDIS,XJACM,XJACI,LNODS,
     .                  POSGP,WEIGP,DERIV,SHAPE,
     .                  COORD,ELCOD,GPCOD)
C***********************************************************************
C
C**** THIS ROUTINE OUTPUTS RESULTS FOR X-Y CURVES
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
      COMMON/TUNING30/PPART,PPARB,PPARI,ESTAB
C
      DIMENSION ELVAR(NSTAT),         eldat(ndata),
     .          DISTO(NTOTV,NDISO),   TLOAD(NTOTV)
      dimension matno(nelem),         proel(nprel,ngrup),
     .          props(nprop,ngrup),   REFOR(NTOTV)
      DIMENSION BUFFE(2),             PRINC(3)
      DIMENSION XJACM(NDIME,*),       XJACI(NDIME,*)
      DIMENSION CARTD(NDIME,NNODE,*), ELDIS(NDOFC,*),
     .          LNODS(NNODE,NELEM)
      DIMENSION POSGP(NDIME,*),       WEIGP(*),
     .          DERIV(NDIME,NNODE,*), SHAPE(NNODE,*)
      DIMENSION COORD(NDIME,*),
     .          ELCOD(NDIME,*),       GPCOD(NDIME,*)
C
      CHARACTER FORM1*13,NLINE*1
      CHARACTER GNUPLO*1
C
      DATA GNUPLO/'#'/
C
      NLINE=CHAR(10)
      FORM1='('//FORMA//',A1)'     ! Format with return
      NRECH=  2                    ! Head record with NPONT
C
      DO ICURV=1,NCURV                                 ! Loop on curves
       NPONT(ICURV+NCOLD)=NPONT(ICURV+NCOLD)+1    ! # of points in curve
       NRECL=5+NPONT(ICURV+NCOLD)                      ! Define record
       DO IAXIS=1,2                                    ! Loop on axes
        ICOMD=MPLOT(ICURV,IAXIS,1)                ! Component & command
        ICOMP=ICOMD/100                           ! Define component
        ICOMD=ICOMD-ICOMP*100                          ! Define command
        IPOI1=MPLOT(ICURV,IAXIS,2)                     ! Define node1
        IPOI2=MPLOT(ICURV,IAXIS,3)                     ! Define node2
        ITOT1=0
        ITOT2=0
        IF(IPOI1.NE.0) ITOT1=(IPOI1-1)*NDOFC+ICOMP     ! Define dof1
        IF(IPOI2.NE.0) ITOT2=(IPOI2-1)*NDOFC+ICOMP     ! Define dof2
        IELEM=IPOI1                                    ! Define element
        IGAUS=IPOI2                                    ! Define Gauss P.
        BVALU=0.0
C
        GOTO(1,2,3,4,5,6,7,8,9,10), ICOMD              ! Do command
C
    1   AVALU=TTIME                                    ! Time
        GO TO 100
    2   AVALU=ISTEP                                    ! Nstep
        GO TO 100
    3   AVALU=TFACT                                    ! Lambda
        GO TO 100
    4   CONTINUE
c       AVALU=TLOAD(ITOT1)
        AVALU=-REFOR(ITOT1)                            ! Reaction
        IF(ITOT2.NE.0)THEN
c        BVALU=TLOAD(ITOT2)
         BVALU=-REFOR(ITOT2)
        ENDIF
        GO TO 100
    5   AVALU=DISTO(ITOT1,1)                           ! Displacement
        IF(ITOT2.NE.0) BVALU=DISTO(ITOT2,1)
        GO TO 100
    6   AVALU=DISTO(ITOT1,2)                           ! Velocity
        IF(ITOT2.NE.0) BVALU=DISTO(ITOT2,2)
        GO TO 100
    7   AVALU=DISTO(ITOT1,3)                           ! Acceleration
        IF(ITOT2.NE.0) BVALU=DISTO(ITOT2,3)
        GO TO 100
    8   CONTINUE
        CALL DATBAS(ELVAR,    3,    2)                 ! Stress
        IF(LARGE.NE.0) THEN
         igrup=matno(ielem)
         NNODL=INT(PROEL(2,IGRUP))
         NNODS=NNODL
         IF(NMEMO5M.EQ.1)
     .    CALL GATHER(DISTO,NDOFC,NPOIN,ELDIS,NDOFC,NNODL,
     .                LNODS(1,IELEM))
         IF(NMEMO1M.EQ.0.OR.NMEMO2M.EQ.0)
     .    CALL DATBAS(ELDAT,    1,    2)          ! geometrical data
         IF(NMEMO1M.EQ.1)
     .    CALL GATHER(COORD,NDIME,NPOIN,ELCOD,NDIME,
     .                NNODL,LNODS(1,IELEM))
         IF(NMEMO2M.EQ.0) THEN
          IF(NMEMO5M.EQ.0) THEN
           CALL PLOCAU(ipoi2,eldat(idata(2)),elvar(istat(1)),
     .                 elvar(istat(4)),XJACM,XJACI,
     .                 POSGP,WEIGP,DERIV,SHAPE,ELCOD,GPCOD)
          ELSE
           CALL PLOCAU(ipoi2,eldat(idata(2)),ELDIS,
     .                 elvar(istat(4)),XJACM,XJACI,
     .                 POSGP,WEIGP,DERIV,SHAPE,ELCOD,GPCOD)
          ENDIF
         ELSE
          CALL PLOCAU(ipoi2,CARTD,ELDIS,
     .                elvar(istat(4)),XJACM,XJACI,
     .                POSGP,WEIGP,DERIV,SHAPE,ELCOD,GPCOD)
         ENDIF                ! nmemo2m.eq.0
        ENDIF                 ! large.ne.0
        igrup=matno(ielem)
        itype=int(proel(5,igrup))
        if(itype.eq.30) then
         imats=int(proel(1,igrup))
         KPARI=INT(PPARI)
         IF(NMEMO1M.EQ.0.OR.NMEMO2M.EQ.0) THEN         ! to be improved
          if(kpari.ne.0) call ssolid(eldat(idata(3)),elvar(istat(4)))
         ENDIF
        endif
c
        ISTAS=ISTAT(4)+(IGAUS-1)*NSTR1
        IF(ICOMP.LE.6) THEN                            !   - cartesian
         AVALU=ELVAR(ISTAS-1+ICOMP)
        ELSE                                           !   - principal
         ICOMP=ICOMP-6
         CALL PRIVAL(NSTR1,ELVAR(ISTAS),PRINC)
         AVALU=PRINC(ICOMP)
        ENDIF
        GO TO 100
    9   CALL DATBAS(ELVAR,    3,    2)                 ! Strain
        ISTAS=ISTAT(3)+(IGAUS-1)*NSTR1
        IF(ICOMP.LE.6) THEN                            !   - cartesian
         AVALU=ELVAR(ISTAS-1+ICOMP)
        ELSE                                           !   - principal
         ICOMP=ICOMP-6
         CALL PRISTN(NSTR1,ELVAR(ISTAS),PRINC)
         AVALU=PRINC(ICOMP)
        ENDIF
        GO TO 100
   10   CALL DATBAS(ELVAR,    3,    2)                 ! Internal Var.
        AVALU=ELVAR(ISTAT(2)+(IGAUS-1)*NHIST-1+ICOMP)
C
  100   BUFFE(IAXIS)=AVALU-BVALU
       ENDDO ! IAXIS
C
       NUNIT=LUCU1-1+ICURV+NCOLD                  ! Define file
       WRITE(UNIT=NUNIT,REC=NRECL,FMT=FORM1)      ! Write record
     .            BUFFE(1),BUFFE(2),NLINE         ! Write head record
C
       IGNUPLO=1                                  ! better as input
       IF(IGNUPLO.EQ.0) THEN
        WRITE(UNIT=NUNIT,REC=NRECH,FMT='(2I5,30X,A1)')
     .                       NPONT(ICURV+NCOLD),0,NLINE
       ELSE
        WRITE(UNIT=NUNIT,REC=NRECH,FMT='(A1,2I5,29X,A1)')
     .                GNUPLO,NPONT(ICURV+NCOLD),0,NLINE
       ENDIF
      ENDDO ! ICURV
C
      RETURN
      END
