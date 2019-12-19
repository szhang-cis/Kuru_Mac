      SUBROUTINE GRAFCARB(TEMP,DELT,FG,FGP,DELFG,XNGN,RGM,FL,FA,FC,MZ,
     .                                 FCP,DELFC,XNCN,RCM,
     .                                 TSOCG,TSOCC,
     .                                 A1,A2,DG,DC,TEG,TEC)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE GRAPHITE AND CARBIDE PHASE-CHANGE
C     FUNCTIONS SIMULTANEOUSLY
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)      
C
C**** GRAPHITE
C
      GU=TEG-TEMP                      ! undercooling
      xngn1=xngn
      delrg=0.0
      IF(TEMP.LT.TEG)THEN
        XNGN1=A1*(GU**2.0)             ! grain density
        DELRG=DG*DELT*(GU**2.0)        ! grain radius increment
C
        if(xngn1.lt.xngn)  then
         xngn1=xngn
        else
         TSOCG=A1*(RGM+DELRG)**3
c        TSOCG=A1*RGM**3               ! other option
        endif
C
        FGP=FG
        FG=12.5664*XNGN1*((RGM+DELRG)**3)/3.0        ! graphite function
        DELFG=FG-FGP
        FL=FL-DELFG
C
        TSOCG=8.0*3.1415/3.0*(TSOCG+XNGN1*(RGM+DELRG)**2*3.0*DG*DELT)*
     .                       (TEG-TEMP)
c       TSOCG=8.0*3.1415/3.0*(TSOCG+XNGN*RGM**2*3.0*DG*DELT)*
c    .                       (TEG-TEMP)              ! other option
      ENDIF
C
C**** CARBIDE
C
      GU=TEC-TEMP      
      xncn1=xncn
      delrc=0.0
      IF(TEMP.LT.TEC)THEN
        XNCN1=A2*(GU**2.0)
        DELRC=DC*DELT*(GU**2.0)
C
        if(xncn1.lt.xncn)  then
         xncn1=xncn
        else
         TSOCC=A2*(RCM+DELRC)**3
c        TSOCC=A2*RCM**3                  ! other option
        endif
C
        FCP=FC
        FC=12.5664*XNCN1*((RCM+DELRC)**3)/3.0
        DELFC=FC-FCP
        FL=FL-DELFC
C
        TSOCC=8.0*3.1415/3.0*(TSOCC+XNCN1*(RCM+DELRC)**2*3.0*DC*DELT)*
     .                       (TEC-TEMP)
c       TSOCC=8.0*3.1415/3.0*(TSOCC+XNCN*RCM**2*3.0*DC*DELT)*
c    .                       (TEC-TEMP)        ! other option
      ENDIF             
C
C**** CORRECTION
C
      IF(FL.LT.0.0)THEN
       alfa=1+fl/(DELFG+DELFC)
       DELFG=DELFG*alfa
       DELFC=DELFC*alfa
       FG=FGP+DELFG
       FC=FCP+DELFC
       FL=0.0
       MZ=100
c
c      the grain density and the grain radius should be corrected
c
      ELSE
       XNGN=XNGN1
       RGM=RGM+DELRG
       XNCN=XNCN1
       RCM=RCM+DELRC
      ENDIF
C
      RETURN
      END
