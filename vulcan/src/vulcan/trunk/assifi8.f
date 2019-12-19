      SUBROUTINE ASSIFI8
C***********************************************************************
C
C**** THIS ROUTINE ASSIGNS MECHANICAL FILES (LINUX)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
C
C***********************************************************************
C
C**** LOGICAL UNITS (LINUX)
C
C***********************************************************************
C
      LUCU1= 191    ! plotting files
      LUCU2= 192
      LUCU3= 193
      LUCU4= 194
      LUCU5= 195
      LUCU6= 196
      LUCU7= 197
      LUCU8= 198
      LUCU9= 199
      LUC10= 200
C
C**** MECHANICAL ASSIGN: INTERNAL FILES
C
      call getenv ('FOR101',ca)
      call getenv ('FOR102',cb) 
      call getenv ('FOR103',cc) 
      call getenv ('FOR104',cd)
      call getenv ('FOR105',ce)
      call getenv ('FOR106',cf)
      call getenv ('FOR107',cg)
      call getenv ('FOR108',ch)
      call getenv ('FOR109',ci)
      call getenv ('FOR110',cj)
      call getenv ('FOR111',ck)
      call getenv ('FOR112',cl)
      call getenv ('FOR113',cm)
      call getenv ('FOR114',cn)
      call getenv ('FOR502',ct)
      call getenv ('FOR139',co)
C
C**** MECHANICAL ASSIGN: EXTERNAL FILES
C
      call getenv ('FOR140',ca1)   ! .geo
      call getenv ('FOR141',cb1)   ! .set
      call getenv ('FOR142',cc1)   ! .mat
      call getenv ('FOR143',cd1)   ! .ini
      call getenv ('FOR144',ce1)   ! .loa
      call getenv ('FOR145',cf1)   ! .fix
      call getenv ('FOR146',cg1)   ! .ini1
      call getenv ('FOR147',ch1)   ! .tun
      call getenv ('FOR148',ci1)   ! .con
      call getenv ('FOR149',cj1)   ! .act
C
C**** MECHANICAL ASSIGN: PLOTTING FILES
C
      call getenv ('FOR191',c1m)
      call getenv ('FOR192',c2m)
      call getenv ('FOR193',c3m)
      call getenv ('FOR194',c4m)
      call getenv ('FOR195',c5m)
      call getenv ('FOR196',c6m)
      call getenv ('FOR197',c7m)
      call getenv ('FOR198',c8m)
      call getenv ('FOR199',c9m)
      call getenv ('FOR200',c10m)
C
      RETURN 
      END
