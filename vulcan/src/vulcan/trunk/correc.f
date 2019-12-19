      SUBROUTINE CORREC(AALPH,DISIT,DISPR,DISTO,DTIME,KINTE,REFOR,
     .                  TALFA,TBETA,TGAMA,TDELT,INDES)
C***********************************************************************
C
C**** THIS ROUTINE UPDATES THE NODAL VARIABLES
C
C     DISIT(1:NTOTV  ) :  nodal iterative   'displacements'
C     DISPR(1:NTOTV,1) :  nodal incremental 'displacements'
C     DISPR(1:NTOTV,2) :  nodal predicted   'velocities'
C     DISPR(1:NTOTV,3) :  nodal predicted   'accelerations'
C     DISTO(1:NTOTV,1) :  nodal current     'displacements'
C     DISTO(1:NTOTV,2) :  nodal current     'velocities'
C     DISTO(1:NTOTV,3) :  nodal current     'accelerations'
C
C     INDES=1 => correction before residual evaluation
C     INDES=2 => correction after residual evaluation
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
C
      DIMENSION DISIT(*), DISPR(NTOTV,*), DISTO(NTOTV,*), REFOR(*)
C
      NPOI1=NPOIN-NPOIC   ! number of points without contact (AL method)
      NTOT1=NPOI1*NDOFC
C
      IF(KDYNA.EQ.1) THEN
       CTIM1=TALFA/(DTIME*DTIME)
       CTIM2=TBETA/DTIME
       CTIM3=TGAMA
      ENDIF
C
      IF(INDES.EQ.1) THEN        ! correction before residual evaluation
C
C**** UPDATE DISPLACEMENTS, VELOCITIES AND ACCELERATIONS
C
       DO ITOTV=1,NTOT1
        REFOR(ITOTV)  =0.0D0
        DISIN         =DISPR(ITOTV,1)+AALPH*DISIT(ITOTV)
C
        IF(KINTE.EQ.0)           ! steady-state
     .   DISTO(ITOTV,1)=DISTO(ITOTV,1)+DISIN
C
        IF(KINTE.EQ.1) THEN      ! Wilson method (at time t+theta dt)
         DISTO(ITOTV,1)=DISTO(ITOTV,1)+DISIN
         DISTO(ITOTV,2)=DISPR(ITOTV,2)+CTIM2*DISIN
         DISTO(ITOTV,3)=DISPR(ITOTV,3)+CTIM1*DISIN
        ENDIF
C
        IF(KINTE.EQ.1.AND.KPROB.EQ.4)
     .   DISTO(ITOTV,3)=DISIN/DTIME      ! Save velocities ??
C
        IF(KINTE.EQ.2) THEN      ! Newmark method
         DISTO(ITOTV,1)=DISTO(ITOTV,1)+DISIN
         DISTO(ITOTV,2)=DISPR(ITOTV,2)+CTIM2*DISIN
         DISTO(ITOTV,3)=DISPR(ITOTV,3)+CTIM1*DISIN
        ENDIF
C
        IF(KINTE.EQ.3) THEN      ! Central difference method
         DISTO(ITOTV,2)=DISPR(ITOTV,2)+CTIM2*DISIN
         DISTO(ITOTV,3)=DISPR(ITOTV,3)+CTIM1*DISIN
        ENDIF
C
        IF(KINTE.EQ.5) THEN      ! HHT method
         DISTO(ITOTV,1)=DISTO(ITOTV,1)+DISIN+(CTIM3-1.0D0)*DISIN
         DISTO(ITOTV,2)=DISPR(ITOTV,2)+CTIM2*DISIN
         DISTO(ITOTV,3)=DISPR(ITOTV,3)+CTIM1*DISIN
        ENDIF
C
        IF(KINTE.EQ.6) THEN      ! Collocation (at time t+theta dt)
         DISTO(ITOTV,1)=DISTO(ITOTV,1)+DISIN
         DISTO(ITOTV,2)=DISPR(ITOTV,2)+CTIM2*DISIN
         DISTO(ITOTV,3)=DISPR(ITOTV,3)+CTIM1*DISIN
        ENDIF
C
       ENDDO
C
       IF(NPOIC.GT.0) THEN
        DO ITOTV=NTOT1+1,NTOTV
         REFOR(ITOTV)  =0.0D0
         DISIN         =DISPR(ITOTV,1)+AALPH*DISIT(ITOTV)
C
C**** INITIALIZATION OF CONTACT FORCE
C
         IF(IAUGM.EQ.0.OR.IAUGM.EQ.1) THEN
          if(nnodc.eq.0) then
cc         IF(DISTO(ITOTV,1).GT.0.0D0) DISTO(ITOTV,1)=0.0D0
ccccc      DISTO(ITOTV,1)=0.0D0          ! to be revised
          else
cc         DISTO(ITOTV,1)=0.0D0          ! to be revised
c          IF(DISTO(ITOTV,1).GT.0.0D0) DISTO(ITOTV,1)=0.0D0
          endif
         ENDIF
         IF(IAUGM.EQ.2.OR.IAUGM.EQ.3) THEN
          IF(DISTO(ITOTV,1).GT.0.0D0) DISTO(ITOTV,1)=0.0D0
         ENDIF
C
         DISTO(ITOTV,1)=DISTO(ITOTV,1)+DISIN
cc       IF(DISTO(ITOTV,1).GT.0.0D0) DISTO(ITOTV,1)=0.0D0
         IF(KINTE.EQ.1) DISTO(ITOTV,2)=DISIN/DTIME
         IF(KINTE.EQ.1.AND.KPROB.EQ.4)
     .    DISTO(ITOTV,3)=DISIN/DTIME     ! Save velocities
         IF(KINTE.EQ.2) THEN
          DISTO(ITOTV,2)=DISPR(ITOTV,2)+CTIM2*DISIN
          DISTO(ITOTV,3)=DISPR(ITOTV,3)+CTIM1*DISIN
         ENDIF
        ENDDO
       ENDIF
      ENDIF
C
      IF(INDES.EQ.2) THEN        ! correction after residual evaluation
       DO ITOTV=1,NTOTV
        DISPR(ITOTV,1)=DISPR(ITOTV,1)+AALPH*DISIT(ITOTV)
C
        IF(KINTE.EQ.1) THEN      ! Wilson method (at time t+dt)
         CONST=3.0D0/TBETA       ! theta
         DISTO(ITOTV,3)=CTIM1/CONST*DISPR(ITOTV,1)-
     .                        CTIM1*DISTO(ITOTV,4)*DTIME+
     .                (1.0D0-TBETA)*DISTO(ITOTV,5)

         DISTO(ITOTV,2)=            DISTO(ITOTV,4)+
     .                        0.5D0*DISTO(ITOTV,3)*DTIME+
     .                        0.5D0*DISTO(ITOTV,5)*DTIME
         DISTO(ITOTV,1)=            DISTO(ITOTV,1)-DISPR(ITOTV,1)+ ! t^U
     .                              DISTO(ITOTV,4)*DTIME+
     .                  1.0D0/6.0D0*DISTO(ITOTV,3)*DTIME*DTIME+
     .                  1.0D0/3.0D0*DISTO(ITOTV,5)*DTIME*DTIME
        ENDIF
C
        IF(KINTE.EQ.3)           ! Central difference method
     .   DISTO(ITOTV,1)=DISTO(ITOTV,1)+DISPR(ITOTV,1)
C
        IF(KINTE.EQ.5)           ! HHT method
     .   DISTO(ITOTV,1)=DISTO(ITOTV,1)-(CTIM3-1.0D0)*DISPR(ITOTV,1)
C
        IF(KINTE.EQ.6) THEN      ! Collocation (at time t+dt)
         CNEW3=TDELT                           ! theta
         CNEW2=1.0D0/(TALFA*CNEW3*CNEW3)       ! alfa
         CNEW1=TBETA*CNEW2*CNEW3               ! delta
         DISTO(ITOTV,3)=CTIM1/CNEW3*DISPR(ITOTV,1)-
     .                        CTIM1*DISTO(ITOTV,4)*DTIME+
     .    (1.0D0-TALFA*CNEW3/2.0D0)*DISTO(ITOTV,5)
         DISTO(ITOTV,2)=            DISTO(ITOTV,4)+
     .                        CNEW1*DISTO(ITOTV,3)*DTIME+
     .                (1.0D0-CNEW1)*DISTO(ITOTV,5)*DTIME
         DISTO(ITOTV,1)=            DISTO(ITOTV,1)-DISPR(ITOTV,1)+ ! t^U
     .                              DISTO(ITOTV,4)*DTIME+
     .                        CNEW2*DISTO(ITOTV,3)*DTIME*DTIME+
     .    (1.0D0-2.0D0*CNEW2)/2.0D0*DISTO(ITOTV,5)*DTIME*DTIME
        ENDIF
C
       END DO
      ENDIF
C
      RETURN
      END
