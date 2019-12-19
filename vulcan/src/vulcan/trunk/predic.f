      SUBROUTINE PREDIC(DISIT,DISPR,DISTO)
C***********************************************************************
C
C**** THIS ROUTINE MAKES A PREDICTION ON NODAL VARIABLES
C
C     DISIT(1:NTOTV  ) :  nodal iterative   'displacements'
C     DISPR(1:NTOTV,1) :  nodal incremental 'displacements'
C     DISPR(1:NTOTV,2) :  nodal predicted   'velocities'
C     DISPR(1:NTOTV,3) :  nodal predicted   'accelerations'
C     DISTO(1:NTOTV,1) :  nodal current     'displacements'
C     DISTO(1:NTOTV,2) :  nodal current     'velocities'
C     DISTO(1:NTOTV,3) :  nodal current     'accelerations'
C
C     DISTO(1:NTOTV,4) :  nodal t-delta t   'displacements'
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION DISIT(*), DISPR(NTOTV,*), DISTO(NTOTV,*)
C
C**** INITIALISE ITERATIVE DISPLACEMENTS
C
      SCALE=0.0D0
      IF(FACPR.NE.0.0D0) SCALE=FACTO/FACPR
C
      IF(LACCE.NE.0.AND.ISTEP.GT.1.AND.SCALE.NE.0.0D0) THEN
       DO ITOTV=1,NTOTV
        DISIT(ITOTV)=DISPR(ITOTV,1)*SCALE
       ENDDO
      ELSE
       DO ITOTV=1,NTOTV
        DISIT(ITOTV)=0.0D0
       ENDDO
      ENDIF
C
C**** INITIALISE INCREMENTAL DISPLACEMENTS  
C
      DO ITOTV=1,NTOTV
       DISPR(ITOTV,1)=0.0D0
C
       IF(KINTE.EQ.1) THEN   ! PREDICT VELOC. & ACCEL.(WILSON'S METHOD)
        DISPR(ITOTV,2)=             -2.0D0*DISTO(ITOTV,2)-
     .                 3.0D0/(TBETA*2.0D0)*DISTO(ITOTV,3)*DTIME
        DISPR(ITOTV,3)=       -2.0D0*TBETA*DISTO(ITOTV,2)/DTIME-
     .                               2.0D0*DISTO(ITOTV,3)
        DISTO(ITOTV,4)=DISTO(ITOTV,2)                    ! t^(dot U)
        DISTO(ITOTV,5)=DISTO(ITOTV,3)                    ! t^(dot dot U)
       ENDIF
C
       IF(KINTE.EQ.2) THEN   ! PREDICT VELOC. & ACCEL.(NEWMARK'S METHOD)
        DISPR(ITOTV,2)=(1.0D0-      TBETA)*DISTO(ITOTV,2)+
     .                 (1.0D0-0.5D0*TBETA)*DISTO(ITOTV,3)*DTIME
        DISPR(ITOTV,3)=            -TALFA *DISTO(ITOTV,2)/DTIME+
     .                 (1.0D0-0.5D0*TALFA)*DISTO(ITOTV,3)
       ENDIF
C
       IF(KINTE.EQ.3) THEN   ! PREDICT VELOC. & ACCEL.(CENT.DIF. METHOD)
        DISPR(ITOTV,2)= TBETA*(DISTO(ITOTV,1)-DISTO(ITOTV,4))/DTIME
        DISPR(ITOTV,3)=-TALFA*(DISTO(ITOTV,1)-DISTO(ITOTV,4))/
     .                                                     (DTIME*DTIME)
        DISTO(ITOTV,4)=DISTO(ITOTV,1)                    ! (t-delta t)^U
       ENDIF
C
       IF(KINTE.EQ.5) THEN   ! PREDICT VELOC. & ACCEL.(HHT'S METHOD)
        DISPR(ITOTV,2)=(1.0D0-      TBETA)*DISTO(ITOTV,2)+
     .                 (1.0D0-0.5D0*TBETA)*DISTO(ITOTV,3)*DTIME
        DISPR(ITOTV,3)=            -TALFA *DISTO(ITOTV,2)/DTIME+
     .                 (1.0D0-0.5D0*TALFA)*DISTO(ITOTV,3)
       ENDIF
C
       IF(KINTE.EQ.6) THEN   ! PREDICT VELOC. & ACCEL.(COLLOCAT. METHOD)
        DISPR(ITOTV,2)=(1.0D0-TBETA*TDELT)*DISTO(ITOTV,2)+
     .     TDELT*(1.0D0-TBETA*TDELT/2.0D0)*DISTO(ITOTV,3)*DTIME
        DISPR(ITOTV,3)=       -TALFA*TDELT*DISTO(ITOTV,2)/DTIME+
     .     (1.0D0-TALFA*TDELT*TDELT/2.0D0)*DISTO(ITOTV,3)
        DISTO(ITOTV,4)=DISTO(ITOTV,2)                    ! t^(dot U)
        DISTO(ITOTV,5)=DISTO(ITOTV,3)                    ! t^(dot dot U)
       ENDIF
C
      ENDDO
C
      RETURN
      END
