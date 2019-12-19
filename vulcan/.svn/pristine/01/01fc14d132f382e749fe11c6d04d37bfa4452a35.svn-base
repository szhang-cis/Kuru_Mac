      SUBROUTINE ASSIFIS8
C***********************************************************************
C
C**** THIS ROUTINE ASSIGNS MICROSTRUCTURAL FILES (LINUX)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
C
C***********************************************************************
C
C**** LOGICAL UNITS (LINUX)
C
C***********************************************************************
C
      LUDATS=315
      LURESS=316
      LUSOLS=317
      LUFRHS=318
      LUPOSS=319
      LUDTSS=320        ! internal files
      LUFROS=321
      LUPRIS=322
      LUSO2S=323
      LUFR2S=324
      LURSTS=325
      LUBFGS=326
      LUPIPS=327
      LUPANS=328
      LUFANS=329

      LUGEOS=350        ! external files
      LUSETS=351
      LUMATS=352
      LUINIS=353
      LULOAS=354
      LUFIXS=355
      LUADVS=356
      LUACTS=359
      LUSTRS=382
C
      LUCU1S=391       ! plotting files
      LUCU2S=392
      LUCU3S=393
      LUCU4S=394
      LUCU5S=395
      LUCU6S=396
      LUCU7S=397
      LUCU8S=398
      LUCU9S=399
      LUC10S=400
C
C**** MICROSTRUCTURAL ASSIGN: INTERNAL FILES
C
      call getenv ('FOR320',cas)
      call getenv ('FOR317',cbs)
      call getenv ('FOR321',ccs)
      call getenv ('FOR318',cds)
      call getenv ('FOR315',ces)
      call getenv ('FOR322',cfs)
      call getenv ('FOR316',cgs)
      call getenv ('FOR323',chs)
      call getenv ('FOR324',cis)
      call getenv ('FOR319',cjs)
      call getenv ('FOR325',cks)
      call getenv ('FOR326',cls)
      call getenv ('FOR327',cms)
      call getenv ('FOR328',cns)
      call getenv ('FOR329',cos)
C
C**** MICROSTRUCTURAL ASSIGN: EXTERNAL FILES
C
      call getenv ('FOR350',ca1s)   ! .geo
      call getenv ('FOR351',cb1s)   ! .set
      call getenv ('FOR352',cc1s)   ! .mat
      call getenv ('FOR353',cd1s)   ! .ini
      call getenv ('FOR354',ce1s)   ! .loa
      call getenv ('FOR355',cf1s)   ! .fix
      call getenv ('FOR356',cg1s)   ! .adv
c
C**** MICROSTRUCTURAL ASSIGN: PLOTTING FILES
C
      call getenv ('FOR391',c1s)
      call getenv ('FOR392',c2s)
      call getenv ('FOR393',c3s)
      call getenv ('FOR394',c4s)
      call getenv ('FOR395',c5s)
      call getenv ('FOR396',c6s)
      call getenv ('FOR397',c7s)
      call getenv ('FOR398',c8s)
      call getenv ('FOR399',c9s)
      call getenv ('FOR400',c10s)
c
C**** COMMON STATUS UNIT
C
      call getenv ('FOR502',cts)
C
      RETURN
      END
