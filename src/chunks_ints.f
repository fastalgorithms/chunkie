      
      subroutine chunkintf(f,fint,ndim,k,nch,chunks,ders,hs)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes the anti-derivative with respect to arc-length 
c     of the real function f (which can be vector valued)
c
c     
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      integer ndim, k, nch
      real *8 f(ndim,k,*), fint(ndim,k,*), chunks(2,k,*), ders(2,k,*)
      real *8 hs(*)
c     local variables
      real *8, allocatable :: ainte(:), adiff(:), x(:), whts(:)
      real *8, allocatable :: ftemp2(:,:), ftemp(:,:), work(:)
      real *8, allocatable :: sums(:), endinter(:)
      integer i, j, l, lwork
      integer korder, itype
      real *8 dsdt

      korder = k-1

      lwork = 3*k**2 + 2*k + 500

      allocate(ainte(k*k),adiff(k*k),x(k),whts(k),ftemp(k,ndim),
     1     ftemp2(k,ndim),work(lwork),sums(ndim),endinter(k))

      itype = 1
      call legeinmt(k,ainte,adiff,x,whts,endinter,itype,work)

      do l = 1,ndim
         sums(l) = 0.0d0
      enddo

      do i = 1,nch
         do j = 1,k
            dsdt = dsqrt(ders(1,j,i)**2+ders(2,j,i)**2)*hs(i)
            do l = 1,ndim
               ftemp2(j,l) = f(l,j,i)*dsdt
            enddo
         enddo
         do l = 1,ndim
c     get spectral anti-derivative
            call amatchunkint(ainte,ftemp2(1,l),ftemp(1,l),k)
         enddo
         do j = 1,k
            do l = 1,ndim
               fint(l,j,i) = ftemp(j,l) + sums(l)
            enddo
         enddo
         do l = 1,ndim
            do j = 1,k
               sums(l) = sums(l) + whts(j)*ftemp2(j,l)
            enddo
         enddo
      enddo

      return
      end

      subroutine amatchunkint(a,x,y,n)
c     evaluates y = a*x
      implicit none
      real *8 a(n,n), x(*), y(*)
      integer n, i, j
      
      do i = 1,n
         y(i) = 0.0d0
      enddo

      do j = 1,n
         do i = 1,n
            y(i) = y(i)+a(i,j)*x(j)
         enddo
      enddo

      return
      end

