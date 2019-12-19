      SUBROUTINE ASSIFIT8
C***********************************************************************
C
C**** THIS ROUTINE ASSIGNS THERMAL FILES (LINUX)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
C
C***********************************************************************
C
C**** LOGICAL UNITS (LINUX)
C
C***********************************************************************
C
      LUDTST=221        ! internal files
      LUSOLT=222
      LUFROT=223
      LUFRHT=224
      LUDATT=225
      LUPRIT=226
      LUREST=227
      LUSO2T=228
      LUFR2T=229
      LUPOST=230
      LURSTT=231
      LUBFGT=232
      LUPIPT=233
      LUPANT=234
      LUFANT=238
C
      LUGEOT=250        ! external files
      LUSETT=251
      LUMATT=252
      LUINIT=253
      LULOAT=254
      LUFIXT=255
      LUADVT=256
      LUACTT=259
      LUSTRT=282
C
      LUCU1T=291       ! plotting files
      LUCU2T=292
      LUCU3T=293
      LUCU4T=294
      LUCU5T=295
      LUCU6T=296
      LUCU7T=297
      LUCU8T=298
      LUCU9T=299
      LUC10T=300
C
C**** THERMAL ASSIGN: INTERNAL FILES
C
      call getenv ('FOR221',cat)
      call getenv ('FOR222',cbt)
      call getenv ('FOR223',cct)
      call getenv ('FOR224',cdt)
      call getenv ('FOR225',cet)
      call getenv ('FOR226',cft)
      call getenv ('FOR227',cgt)
      call getenv ('FOR228',cht)
      call getenv ('FOR229',cit)
      call getenv ('FOR230',cjt)
      call getenv ('FOR231',ckt)
      call getenv ('FOR232',clt)
      call getenv ('FOR233',cmt)
      call getenv ('FOR234',cnt)
      call getenv ('FOR238',cot)
C
C**** THERMAL ASSIGN: EXTERNAL FILES
C
      call getenv ('FOR250',ca1t)   ! .geo
      call getenv ('FOR251',cb1t)   ! .set
      call getenv ('FOR252',cc1t)   ! .mat
      call getenv ('FOR253',cd1t)   ! .ini
      call getenv ('FOR254',ce1t)   ! .loa
      call getenv ('FOR255',cf1t)   ! .fix
      call getenv ('FOR256',cg1t)   ! .adv
      call getenv ('FOR259',ch1t)   ! .act
      call getenv ('FOR282',ci1t)   ! .str
C
C**** THERMAL ASSIGN: PLOTTING FILES
C
      call getenv ('FOR291',c1t)
      call getenv ('FOR292',c2t)
      call getenv ('FOR293',c3t)
      call getenv ('FOR294',c4t)
      call getenv ('FOR295',c5t)
      call getenv ('FOR296',c6t)
      call getenv ('FOR297',c7t)
      call getenv ('FOR298',c8t)
      call getenv ('FOR299',c9t)
      call getenv ('FOR300',c10t)
c
C**** COMMON STATUS UNIT
C
      call getenv ('FOR502',ctt)
C
      RETURN
      END
