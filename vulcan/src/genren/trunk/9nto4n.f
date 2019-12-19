      program c9nto4n          
C***********************************************************************
C
C PROGRAMA LECTOR DE ELEMENTOS DE NUEVE NODOS
C LOS TRANSFORMA EN ELEMENTOS DE CUATRO NODOS
C 
C***********************************************************************
      implicit double precision(a-h,o-z)
      character*11 filein
      character*11 fileout

      dimension nods(9)

      write(6,*)'NOMBRE DEL INPUT (11 ch.)='
      read(5,*)filein
      write(6,*)'NOMBRE DEL OUTPUT ='
      read(5,*)fileout

      open (3,file=filein)
      open (4,file=fileout)

      kelem=0

      read(3,201)nelem
201   format(i5)

      do i=1,nelem

      read(3,202)ielem,(nods(l),l=1,9)                       
202   format(10i5)

      kelem=kelem+1 
      write(4,203)kelem,nods(1),nods(5),nods(9),nods(8)
      kelem=kelem+1 
      write(4,203)kelem,nods(5),nods(2),nods(6),nods(9)
      kelem=kelem+1 
      write(4,203)kelem,nods(8),nods(9),nods(7),nods(4)
      kelem=kelem+1 
      write(4,203)kelem,nods(9),nods(6),nods(3),nods(7)
203   format(i5,4(1x,i4))

      end do

      end
