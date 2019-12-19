      SUBROUTINE ASSIFIS3
C***********************************************************************
C
C**** THIS ROUTINE ASSIGNS MICROSTRUCTURAL FILES (VAX)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'

      call runends('error: vax not implemented')

C
C**** COUPLING ASSIGN              ! it should be in assifid !!!!!
C
cts      call getenv ('FOR037',ccoc)
cts      call getenv ('FOR237',ccod) 
C
C**** MICROSTRUCTURAL ASSIGN: INTERNAL FILES
C
cts      call getenv ('FOR021',cas)
cts      call getenv ('FOR022',cbs) 
cts      call getenv ('FOR023',ccs) 
cts      call getenv ('FOR024',cds)
          call getenv ('FOR025',ces)
cts      call getenv ('FOR026',cfs)
          call getenv ('FOR027',cgs)
cts      call getenv ('FOR028',chs)
cts      call getenv ('FOR029',cis)
cts      call getenv ('FOR030',cjs)
cts      call getenv ('FOR031',cks)
cts      call getenv ('FOR032',cls)
cts      call getenv ('FOR033',cms)
C
C**** MICROSTRUCTURAL ASSIGN: EXTERNAL FILES
C
cts      call getenv ('FOR050',ca1s)   ! .geo
cts      call getenv ('FOR051',cb1s)   ! .set
cts      call getenv ('FOR052',cc1s)   ! .mat
cts      call getenv ('FOR053',cd1s)   ! .ini
cts      call getenv ('FOR054',ce1s)   ! .loa
cts      call getenv ('FOR055',cf1s)   ! .fix
C
C**** MICROSTRUCTURAL ASSIGN: PLOTTING FILES
C
cts      call getenv ('FOR121',c1ms)
cts      call getenv ('FOR122',c2ms)
cts      call getenv ('FOR123',c3ms)
cts      call getenv ('FOR124',c4ms)
cts      call getenv ('FOR125',c5ms)
cts      call getenv ('FOR126',c6ms)
cts      call getenv ('FOR127',c7ms)
cts      call getenv ('FOR128',c8ms)
cts      call getenv ('FOR129',c9ms)
cts      call getenv ('FOR130',c10ms)
C
      RETURN
      END
