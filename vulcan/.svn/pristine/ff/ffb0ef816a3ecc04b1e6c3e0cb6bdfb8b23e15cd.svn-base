c=======================================================================
      subroutine gmrespt(a,diag, idim,work, x, k0, 
     .                         itmax, toler, yours, iprint,
     .                         UNKNW,RHSID,LNUEQ,LNODS,NDOFN,NELEM,
     .                         NEVAB,NNODE,IPASS,NTOTV,NPOIN,
     .                         LURES)
c-----------------------------------------------------------------------
c     Implemented by: Didier Joannas
c                     (modified by j.c.: preconditioned GMRES)  
c                     (modified by YO)
c    DATAS :
c            idim   : dimension of the problem (NEQNS)
c            k0     : dimension of the krylov subspace
c            itmax  : maximal number of iterations
c            toler  : relative precision needed to stop the process
c            Ax     : external procedure wich computes Ax or b-Ax
c            iprint : = 0 nothing is printed
c                     = 1 the initial residual and the final relative
c                         residual are printed
c                     = 2 the inial residual and all the relative
c                         residuals are printed
c                                       
c    RESULT :
c            x : solution of Ax=b
c                                                
c    WORK ARRAYS :
c            work  : it's a work array I need for gmres
c                    the dimension of this array must be, at least of
c                    idim*(k0+1)+(k0*(k0+1))/2+4*k0+2
c            yours : it's a work array you may need to compute Ax or
c                    b-Ax
c                                                   
c     description of the Ax subroutine :
c         subroutine Ax(x,idim,result,yours,job)
c         job = 0 : compute Ax,
c             = 1 : compute b-Ax
c         the result is put in the vector result
c
c     VERSION  2.0   Oct - 93     
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)        
      integer*4 idim, k0, itmax, iprint
      integer*4 lu, luj, luipl1, lsin, lcos, le, ly, lwork
      integer*4 lh, lhi, lhip1i, lhji, iter
      integer*4 iincs
      dimension x(*), work(*), yours(*)
      dimension diag(*)
C
      DIMENSION UNKNW(*),RHSID(*),LNUEQ(*),LNODS(*)
      INTEGER*4 NDOFN,NELEM,NEVAB,NNODE,IPASS,NTOTV,NPOIN,
     .          LURES
C

c     write(lures,*) '==========================================='
c     if (iprint.gt.0) then                           ! original version
       write(lures,*)'Welcome in gmres subroutine....'
       write(lures,*)'equations',idim,'  Krylov dim:',k0
c     endif

c     the first approximation of the solution is zero :
      do  i = 1, idim
       x(i)=0.d0
      end do

c     we reserve place for the basis of krylov subspaces, for the
c     hessenberg matrix, for cos(teta) and sin(teta), for y and e
c     matrix, for cos(teta) and sin(teta), for y and e
      lu   = 1
      lh   = lu+idim*(k0+1)
      lcos = lh+(k0*(k0+1))/2+1
      lsin = lcos+k0
      le   = lsin+k0
      ly   = le+k0+1
      lwork= ly+k0
      iter = 0
      KPRIM=1
c     external loop
 20   continue
c     compute the first residual from zero or from the last solution got
      call GMRESAT(a,x,diag,work(lu),yours,NPOIN,1,LNUEQ,LNODS,UNKNW,
     .             RHSID,NDOFN,NELEM,NEVAB,NNODE,IPASS,NTOTV,KPRIM,
     .             idim)                             ! added variable
                                                     ! preconditioner...
      call scalpr(work(lu),work(lu),idim,prod)
      prod   = dsqrt(prod)
      prodin = 1.d0/prod
      call sctvec(prodin,work(lu),work(lu),idim)
      work(le) = prod
c     if ((iprint.gt.0).and.(iter.eq.0)) then         ! original version
      if (                  (iter.eq.0)) then
       write(lures,*)'initial residual :              ',prod
       rinit = prod
      endif

c     initialisation of the error vector :
      do i = 1, k0
       work(le+i) = 0.d0
      end do

      lhi = lh
      i = 0

c     internal loop : minimisation of the residual on a krylov subspace
35    i = i+1
      luipl1 = lu+i*idim 
      luipl2=luipl1-idim
      call GMRESAT(a,work(luipl2),diag,work(luipl1),
     .             yours,NPOIN,0,LNUEQ,LNODS,UNKNW,
     .             RHSID,NDOFN,NELEM,NEVAB,NNODE,IPASS,NTOTV,KPRIM,
     .             idim)                             ! added variable
                                                     ! preconditioner...
c     do iii = 1 , idim  
c      write(lures,*) work(luipl1-idim+ iii-1),work(luipl1 + iii-1),
c    .            yours(iii),' d',diag(iii)
c     end do 
c     write(lures,*)  '=================='
c   gram-schmidt orthogonalisation and construction of the column number
c   i of the hessenberg matrix :
      do j = 1, i
       luj  = lu+(j-1)*idim
       lhji = lhi+j-1
       call scalpr(work(luipl1),work(luj),idim,work(lhji))
       coef = -work(lhji)
       call advect(1.d0,work(luipl1),coef,work(luj),work(luipl1),idim)
      end do

      lhip1i = lhi+i  
      call scalpr(work(luipl1),work(luipl1),idim,work(lhip1i))
      work(lhip1i) = dsqrt(work(lhip1i))
      
      coef = 1.d0/work(lhip1i)
      call sctvec(coef,work(luipl1),work(luipl1),idim)
                                     
c     QR factorisation of the hessenberg matrix                      
      do j = 1, i-1
       lhji         = lhi+j-1
       lhjp1i       = lhi+j
       hji          = work(lhji)   
       hjp1i        = work(lhjp1i) 
       work(lhji)   = work(lcos+j-1)*hji+work(lsin+j-1)*hjp1i
       work(lhjp1i) = work(lcos+j-1)*hjp1i-work(lsin+j-1)*hji
      end do          
      hii            = work(lhi+i-1)
      hip1i          = work(lhi+i)  
      sqrhi          = dsqrt(hii*hii+hip1i*hip1i)
      work(lcos+i-1) = hii/sqrhi
      work(lsin+i-1) = hip1i/sqrhi
      work(lhi+i-1)  = sqrhi
cc                            work(lhi+i)    = 0.d0
      ei             = work(le+i-1)
      work(le+i-1)   = work(lcos+i-1)*ei
      work(le+i)     = -work(lsin+i-1)*ei
      residu         = abs(work(le+i)/rinit)
      iter = iter+1
c     write(lures,*)                                  ! original version
c    . 'relative residual number',iter,'=',residu
      if((iprint.gt.0).or.(iprint.eq.0.and.iter.eq.1))
     . write(lures,*)'relative residual number',iter,'=',residu
      lhi = lhi+i
      if (((residu.gt.toler).and.(iter.lt.itmax)).and.(i.lt.k0)) goto 35
c                 end of the internal loop
           
c     resolution of the upper triangular system :
      call solve(work(lh),work(ly),work(le),i)
                        
c     computation of the new estimation of the solution :
      do j = 1, i
       luj= lu+(j-1)*idim
       call advect(1.d0,x,work(ly+j-1),work(luj),x,idim)
      end do
                              
c     test for reinitialization :
      if ((residu.gt.toler).and.(iter.lt.itmax)) goto 20
c     if (iprint.gt.0)                                ! original version
c    . write(lures,*)'final relative residual =',residu
      write(lures,*)'final relative residual',iter,'=',residu
      if ((residu.gt.toler).and.(iter.ge.itmax))
     .   write(lures,*) 'salio por maximo numero de iteraciones' 
      if ((residu.gt.toler).and.(iter.ge.itmax)) then
         write(lures,*)
     .        'Exited of GMRES solver. Maximum iteration number reached'
         write(lures,*)
     .        'Last rel. error obtained was:',residu,' iter:',iter
         write(lures,*)
     .   'Maximum iteration number.. or.. 
     .    krylov dimension should be increased'

      end if
      return
      end
