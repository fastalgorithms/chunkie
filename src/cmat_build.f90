
        !
        !     This file contains several subroutines for constructing
        !     matrices corresponding to layer potentials with
        !     complex-valued kernels, i.e. helmholtz, etc.
        !
        !     The important subroutines are as follows:
        !
        !       zbuildmat - 
        !
        !

      subroutine zget_neumann(wgeo, fincoming, par0, src, pars1, &
          pars2, uinc)
        implicit real *8 (a-h,o-z)
        real *8, intent(in) :: wgeo(*), src(2)
        complex *16, intent(in) :: par0, pars1(*), pars2(*)
        complex *16, intent(out) :: uinc(*)

        real *8 :: dnorm(10)
        complex *16 :: val, grad(2), hess(2,2)

        !
        ! evaluate the dirichlet data on the boundary contained
        ! in wgeo using the function fincoming
        !

        call chunkunpack1(wgeo, k, nch, ichunks, iadjs, iders, &
            iders2, ihs)

        do i = 1,k*nch
          ind = 2*(i-1)
          xp = wgeo(iders+ind)
          yp = wgeo(iders+ind+1)
          sc = sqrt(xp**2+yp**2)
          dnorm(1) = yp/sc
          dnorm(2) = -xp/sc
          call fincoming(par0, src, wgeo(ichunks+ind), pars1, &
              pars2, val, grad, hess)
          uinc(i) = (dnorm(1)*grad(1) + dnorm(2)*grad(2))
        enddo

      end subroutine




      subroutine zget_dirichlet(wgeo, fincoming, par0, src, pars1, &
          pars2, uinc)
        implicit real *8 (a-h,o-z)
        real *8, intent(in) :: wgeo(*), src(2)
        complex *16, intent(in) :: par0, pars1(*), pars2(*)
        complex *16, intent(out) :: uinc(*)

        complex *16 :: val, grad(2), hess(2,2)

        !
        ! evaluate the dirichlet data on the boundary contained
        ! in wgeo using the function fincoming
        !

        call chunkunpack1(wgeo, k, nch, ichunks, iadjs, iders, &
            iders2, ihs)

        do i = 1,k*nch
          ind = 2*(i-1)
          call fincoming(par0, src, wgeo(ichunks+ind), pars1, &
              pars2, val, grad, hess)
          uinc(i) = val
        enddo

      end subroutine





      subroutine zkernel_cfie(par0, src, targ, src_normal, &
        targ_normal, alpha, beta, fgreen, pars1, pars2, val)
        implicit real *8 (a-h,o-z)
        real *8, intent(in) :: src(2), targ(2), src_normal(2)
        real *8, intent(in) :: targ_normal(2)
        complex *16, intent(in) :: par0, alpha, beta, pars1(*)
        complex *16, intent(in) :: pars2(*)
        complex *16, intent(out) :: val

        complex *16, parameter :: ima = (0,1)
        complex *16 :: val2, grad2(2), hess2(2,2), dlp, slp
        !
        ! this computes the combined field equation for the exterior
        ! dirichlet problem, namely:
        !
        !         alpha*S + beta*D
        !

        call fgreen(par0, src, targ, pars1, pars2, val2, &
            grad2, hess2)
        dlp = -(src_normal(1)*grad2(1) + src_normal(2)*grad2(2))
        slp = val2
        val = alpha*slp + beta*dlp

        return
      end subroutine





      subroutine zkernel_cfie_hyper(par0, src, targ, src_normal, &
        targ_normal, alpha, beta, fgreen, pars1, pars2, val)
        implicit real *8 (a-h,o-z)
        real *8, intent(in) :: src(2), targ(2), src_normal(2)
        real *8, intent(in) :: targ_normal(2)
        complex *16, intent(in) :: par0, alpha, beta, pars1(*)
        complex *16, intent(in) :: pars2(*)
        complex *16, intent(out) :: val

        complex *16, parameter :: ima = (0,1)
        complex *16 :: val2, grad2(2), hess2(2,2), dlp, slp, vec(10)
        !
        ! this computes the combined field equation for (what is
        ! usually) the transmission problem - note that unless this
        ! routine receives a difference kernel, the operators are
        ! probably quite singular. The matrix build is:
        !
        !         alpha*S' + beta*D'
        !
        ! where ' denotes normal derivative at the target
        !

        call fgreen(par0, src, targ, pars1, pars2, val2, &
            grad2, hess2)

        slp = grad2(1)*targ_normal(1) + grad2(2)*targ_normal(2)       

        vec(1) = hess2(1,1)*src_normal(1) + hess2(1,2)*src_normal(2)
        vec(2) = hess2(2,1)*src_normal(1) + hess2(2,2)*src_normal(2)
        dlp = -(vec(1)*targ_normal(1) + vec(2)*targ_normal(2))

        val = alpha*slp + beta*dlp

        return
      end subroutine





      subroutine zkernel_slp(par0, src, targ, src_normal, &
          targ_normal, q1, q2, fgreen, pars1, pars2, val)
        implicit real *8 (a-h,o-z)
        real *8, intent(in) :: src(2), targ(2), src_normal(2)
        real *8, intent(in) :: targ_normal(2), pars1(*), pars2(*)
        complex *16, intent(in) :: par0, q1, q2
        complex *16, intent(out) :: val

        complex *16 :: val2, grad2(2), hess2(2,2)
        !
        ! this is just a wrapper routine for the double layer
        ! kernel evaluation using fgreen
        !

        call fgreen(par0, src, targ, pars1, pars2, val, &
            grad2, hess2)

        return
      end subroutine





      subroutine zkernel_dlp(par0, src, targ, src_normal, &
          targ_normal, q1, q2, fgreen, pars1, pars2, val)
        implicit real *8 (a-h,o-z)
        real *8, intent(in) :: src(2), targ(2), src_normal(2)
        real *8, intent(in) :: targ_normal(2), pars1(*), pars2(*)
        complex *16, intent(in) :: par0, q1, q2
        complex *16, intent(out) :: val

        complex *16 :: val2, grad2(2), hess2(2,2)
        !
        ! this is just a wrapper routine for the double layer
        ! kernel evaluation using fgreen
        !

        call fgreen(par0, src, targ, pars1, pars2, val2, &
            grad2, hess2)
        val = -(src_normal(1)*grad2(1) + src_normal(2)*grad2(2))

        return
      end subroutine





      subroutine zkernel_dprime(par0, src, targ, src_normal, &
          targ_normal, q1, q2, fgreen, pars1, pars2, val)
        implicit real *8 (a-h,o-z)
        real *8, intent(in) :: src(2), targ(2), src_normal(2)
        real *8, intent(in) :: targ_normal(2), pars1(*), pars2(*)
        complex *16, intent(in) :: par0, q1, q2
        complex *16, intent(out) :: val

        complex *16 :: val2, grad2(10), hess2(2,2), vec(10)
        !
        ! this is just a wrapper routine for the kernel of the 
        ! normal derivative of the double layer - which is a hyper
        ! singular operator, unless the kernel is a difference kernel
        !

        call fgreen(par0, src, targ, pars1, pars2, val2, &
            grad2, hess2)

        vec(1) = hess2(1,1)*src_normal(1) + hess2(1,2)*src_normal(2)
        vec(2) = hess2(2,1)*src_normal(1) + hess2(2,2)*src_normal(2)

        val = -(vec(1)*targ_normal(1) + vec(2)*targ_normal(2))
        
        return
      end subroutine







      subroutine zkernel_sprime(par0, src, targ, src_normal, &
          targ_normal, q1, q2, fgreen, pars1, pars2, val)
        implicit real *8 (a-h,o-z)
        real *8, intent(in) :: src(2), targ(2), src_normal(2)
        real *8, intent(in) :: targ_normal(2), pars1(*), pars2(*)
        complex *16, intent(in) :: par0, q1, q2
        complex *16, intent(out) :: val

        complex *16 :: val2, grad2(2), hess2(2,2)
        !
        ! this is just a wrapper routine for the kernel of the 
        ! normal derivative of the single layer
        !

        call fgreen(par0, src, targ, pars1, pars2, val2, &
            grad2, hess2)
        val = targ_normal(1)*grad2(1) + targ_normal(2)*grad2(2)

        return
      end subroutine

      subroutine chunks_dcopy(n, x, y)
        implicit real *8 (a-h,o-z)
        real *8 :: x(n), y(n)

        do i = 1,n
           y(i) = x(i)
        enddo

        return
      end subroutine chunks_dcopy





      subroutine zmat_by_dmat(m, n, a, k, b, c)
        implicit real *8 (a-h,o-z)
        real *8 :: b(n,k)
        complex *16 :: a(m,n), c(m,k), dd

        !
        ! this routine multiplies a complex matrix by a real matrix
        !
        
        do i = 1,m
           do j = 1,k
              dd = 0
              do ijk = 1,n
                 dd = dd + a(i,ijk)*b(ijk,j)
              enddo
              c(i,j) = dd
           enddo
        enddo

        return
      end subroutine





      subroutine dmat_by_zmat(m, n, a, k, b, c)
        implicit real *8 (a-h,o-z)
        real *8 :: a(m,n)
        complex *16 :: b(n,k), c(m,k), dd

        !
        ! this routine multiplies a complex matrix by a real matrix
        !
        
        do i = 1,m
           do j = 1,k
              dd = 0
              do ijk = 1,n
                 dd = dd + a(i,ijk)*b(ijk,j)
              enddo
              c(i,j) = dd
           enddo
        enddo

        return
      end subroutine


      subroutine blas2(m, n, a, x, y)
        implicit real *8 (a-h,o-z)
        real *8 :: a(m,n), x(n), y(m)
        character *1 :: trans

        trans = 'n'
        alpha = 1
        lda = m
        incx = 1
        beta = 0
        incy = 1

        call dgemv(trans, m, n, alpha, a, lda, x, &
            incx, beta, y, incy)

        return
      end subroutine





      subroutine zblas2(m, n, a, x, y)
        implicit real *8 (a-h,o-z)
        complex *16 :: a(m,n), x(n), y(m)
        character *1 :: trans

        trans = 'n'
        alpha = 1
        lda = m
        incx = 1
        beta = 0
        incy = 1

        call zgemv(trans, m, n, alpha, a, lda, x, &
            incx, beta, y, incy)

        return
      end subroutine





      subroutine zapplymat(norder, wgeo, fkernel, q1, q2, &
          fgreens, par0, pars1, pars2, ifself, ifnear, iffar, &
          sigma, pot)
        implicit real *8 (a-h,o-z)
        real *8 :: wgeo(*)
        complex *16 :: q1(*), q2(*), par0, pars1(*), pars2(*)
        complex *16 :: sigma(*), pot(*)
        external fkernel, fgreens

        !
        ! this routine is analogous to zbuildmat, but only applies
        ! the operator (in blocks)
        !

        call chunkunpack1(wgeo, k, nch, ichunks, iadjs, iders, &
            iders2, ihs)

        ntot = k*nch
        if (norder .ne. k) then
          call prinf('norder = *', norder, 1)
          call prinf('k = *', k, 1)
          call prinf('bomb!!!!*', k, 0)
          stop
        endif

        call zapplymat0(norder, k, nch, wgeo, fkernel, q1, q2, &
            fgreens, par0, pars1, pars2, ifself, ifnear, iffar, &
            sigma, pot)

        return
      end subroutine


      subroutine zapplymat0(norder, k, nch, wgeo, fkernel, q1, q2, &
          fgreens, par0, pars1, pars2, ifself, ifnear, iffar, &
          sigma, pot)
        real *8 :: wgeo(*)
        complex *16 :: q1(*), q2(*), par0, pars1(*), pars2(*)
        complex *16 :: sigma(k,nch), pot(k,nch)

        external fkernel, fgreens

        real *8, allocatable :: xs(:), whts(:)
        real *8, allocatable :: xs1(:), whts1(:)
        real *8, allocatable :: xs0(:,:), whts0(:,:)
        complex *16 :: pottemp(10000)
        complex *16, allocatable :: tempmat(:,:)

        !
        ! apply the operator fkernel(fgreens) piece by piece
        !
        done = 1
        allocate(tempmat(k,k))

        itype = 1
        allocate(xs(k), whts(k))
        call legerts(itype, k, xs, whts)

        !
        ! get the bremer near and self quadratures
        !
        call getquadsinfo(k, nquad1, nquad0)
        allocate(xs1(nquad1), whts1(nquad1))
        allocate(xs0(k,k), whts0(k,k))
        call getquads(k, nquad1, xs1, whts1, nquad0, xs0, whts0)

        !
        ! initialize the potential
        !
        do i = 1,nch
          do j =1,k
            pot(j,i) = 0
          enddo
        enddo

        do itarget = 1,nch
          imat = 1 + (itarget-1)*k

          call chunk_neighbors(itarget, wgeo, ibefore, iafter)

          !
          ! for all target chunks, loop over source chunks
          !
          do jsource = 1,nch
            jmat = 1 + (jsource-1)*k

            if ((jsource .eq. ibefore) .and. (ifnear .ne. 0)) then
              call znear_buildmat(itarget, jsource, k, wgeo, &
                  fkernel, q1, q2, fgreens, par0, pars1, pars2, &
                  nquad1, xs1, whts1, tempmat)
              call zblas2(k, k, tempmat, sigma(1,jsource), pottemp)
              call addpot(k, pottemp, pot(1,itarget))
              cycle
            endif

            if ((jsource .eq. iafter) .and. (ifnear .ne. 0)) then
              call znear_buildmat(itarget, jsource, k, wgeo, &
                  fkernel, q1, q2, fgreens, par0, pars1, pars2, &
                  nquad1, xs1, whts1, tempmat)
              call zblas2(k, k, tempmat, sigma(1,jsource), pottemp)
              call addpot(k, pottemp, pot(1,itarget))
              cycle
            endif

            if ((jsource .eq. itarget) .and. (ifself .ne. 0)) then
              call zdiagonal_buildmat(itarget, k, wgeo, fkernel, &
                  q1, q2, fgreens, par0, pars1, pars2, xs0, &
                  whts0, tempmat)
              call zblas2(k, k, tempmat, sigma(1,jsource), pottemp)
              call addpot(k, pottemp, pot(1,itarget))
              cycle
            endif

            if (iffar .ne. 0) then
              call zfar_buildmat(itarget, jsource, k, wgeo, &
                  fkernel, q1, q2, fgreens, par0, pars1, pars2, xs, &
                  whts, tempmat)
              call zblas2(k, k, tempmat, sigma(1,jsource), pottemp)
              call addpot(k, pottemp, pot(1,itarget))
            endif
            
          enddo
        enddo

      end subroutine





      subroutine zapply_self(norder, k, nch, ich, wgeo, fkernel, &
          q1, q2, fgreens, par0, pars1, pars2, sigma, pot)
        real *8 :: wgeo(*)
        complex *16 :: q1(*), q2(*), par0, pars1(*), pars2(*)
        complex *16 :: sigma(k,nch), pot(k,nch)

        external fkernel, fgreens

        real *8, allocatable :: xs(:), whts(:)
        real *8, allocatable :: xs1(:), whts1(:)
        real *8, allocatable :: xs0(:,:), whts0(:,:)
        complex *16, allocatable :: tempmat(:,:)

        !
        ! apply the operator fkernel(fgreens) piece by piece
        !
        done = 1
        allocate(tempmat(k,k))

        itype = 1
        allocate(xs(k), whts(k))
        call legerts(itype, k, xs, whts)

        !
        ! get the bremer near and self quadratures
        !
        call getquadsinfo(k, nquad1, nquad0)
        allocate(xs1(nquad1), whts1(nquad1))
        allocate(xs0(k,k), whts0(k,k))
        call getquads(k, nquad1, xs1, whts1, nquad0, xs0, whts0)

        call zdiagonal_buildmat(ich, k, wgeo, fkernel, &
            q1, q2, fgreens, par0, pars1, pars2, xs0, &
            whts0, tempmat)
        call zblas2(k, k, tempmat, sigma(1,ich), pot(1,ich))

        return
      end subroutine




      subroutine zapply_near(norder, k, nch, ich, ich_src, wgeo, &
          fkernel, q1, q2, fgreens, par0, pars1, pars2, sigma, pot)
        real *8 :: wgeo(*)
        complex *16 :: q1(*), q2(*), par0, pars1(*), pars2(*)
        complex *16 :: sigma(k,nch), pot(k,nch)

        external fkernel, fgreens

        real *8, allocatable :: xs(:), whts(:)
        real *8, allocatable :: xs1(:), whts1(:)
        real *8, allocatable :: xs0(:,:), whts0(:,:)
        complex *16, allocatable :: tempmat(:,:)

        !
        ! apply the operator fkernel(fgreens) piece by piece
        !
        done = 1
        allocate(tempmat(k,k))

        itype = 1
        allocate(xs(k), whts(k))
        call legerts(itype, k, xs, whts)

        !
        ! get the bremer near and self quadratures
        !
        call getquadsinfo(k, nquad1, nquad0)
        allocate(xs1(nquad1), whts1(nquad1))
        allocate(xs0(k,k), whts0(k,k))
        call getquads(k, nquad1, xs1, whts1, nquad0, xs0, whts0)

        call znear_buildmat(ich, ich_src, k, wgeo, &
            fkernel, q1, q2, fgreens, par0, pars1, pars2, &
            nquad1, xs1, whts1, tempmat)
        call zblas2(k, k, tempmat, sigma(1,ich_src), pot(1,ich))

        return
      end subroutine

      subroutine addpot(k, pot1, pot2)
        implicit real *8 (a-h,o-z)
        complex *16 :: pot1(k), pot2(k)

        do i = 1,k
          pot2(i) = pot2(i) + pot1(i)
        enddo

        return
      end subroutine

      
      subroutine zbuildmat(norder, wgeo, fkernel, q1, q2, &
          fgreens, par0, pars1, pars2, ntot, cmat)
        implicit real *8 (a-h,o-z)
        real *8, intent(in) :: wgeo(*)
        complex *16, intent(in) :: par0, pars1(*), pars2(*)
        complex *16, intent(in) :: q1(*), q2(*)
        integer, intent(out) :: ntot
        complex *16, intent(out) :: cmat(*)
        external fkernel, fgreens

        !
        ! build **an arbitrary** layer potential matrix using bremer-style
        ! chunk quadratures. this assumed that each chunk is *locally*
        ! oriented counter-clockwise. 
        !
        ! input:
        !   norder - the order of the quadrature routine to use even
        !       though right now it enforces that norder is equal to
        !       the number of points on each chunk
        !   wgeo - the geometry structure created by chunks routines
        !   fkernel - the kernel evaluation routine, must have its
        !       calling sequence of the form:
        !           fkernel(par0, src, targ, src_normal, targ_normal,
        !               q1, q2, fgreens, pars1, pars2, val)
        !   fgreens- subroutine generating the greens function, its
        !       calling sequence must be of the form:
        !           fgreens(par0, src, targ, pars1, pars2, val, grad,
        !               hess)
        !     NOTE: be sure to take note that fgreens is **passed** to
        !         to the subroutine fkernel - that is, we make a
        !         distinction between the greens function and the kernel
        !
        ! output:
        !   ntot - the dimension of the matrix, ntot = k*nch
        !   cmat - the matrix which applies the kernel fkernel
        !
        !     This is mainly a memory management and unpacking routine
        !

        !
        ! first get some chunk info
        !
        call chunkunpack1(wgeo, k, nch, ichunks, iadjs, iders, &
            iders2, ihs)

        ntot = k*nch
        if (norder .ne. k) then
          call prinf('norder = *', norder, 1)
          call prinf('k = *', k, 1)
          call prinf('bomb!!!!*', k, 0)
          stop
        endif

        call zbuildmat0(norder, k, nch, wgeo, fkernel, q1, q2, &
            fgreens, par0, pars1, pars2, ntot, cmat)

        return
      end subroutine





      subroutine zbuildmat0(norder, k, nch, wgeo, fkernel, q1, q2, &
          fgreens, par0, pars1, pars2, ntot, zmat)
        implicit real *8 (a-h,o-z)
        integer, intent(in) :: norder, k, nch, ntot
        real *8, intent(in) :: wgeo(*)
        complex *16, intent(in) :: par0, pars1(*), pars2(*)
        complex *16, intent(in) :: q1(*), q2(*)
        complex *16, intent(out) :: zmat(ntot,ntot)

        external fkernel, fgreens

        real *8, allocatable :: xs(:), whts(:)
        real *8, allocatable :: xs1(:), whts1(:)
        real *8, allocatable :: xs0(:,:), whts0(:,:)
        !complex *16, allocatable :: tempmat(:,:)
        complex *16 :: tempmat(10000)

        !
        ! actually go ahead and build the matrix, block by block
        !
        done = 1
        !allocate(tempmat(k,k))

        itype = 1
        allocate(xs(k), whts(k))
        call legerts(itype, k, xs, whts)

        !
        ! get the bremer near and self quadratures
        !

        do j = 1,ntot
          do i = 1,ntot
            zmat(i,j) = 0
          enddo
        enddo
            
        call getquadsinfo(k, nquad1, nquad0)
        allocate(xs1(nquad1), whts1(nquad1))
        allocate(xs0(k,k), whts0(k,k))
        call getquads(k, nquad1, xs1, whts1, nquad0, xs0, whts0)

        !$omp parallel do default(shared) &
        !$omp     private(itarget, imat, ibefore, iafter) &
        !$omp     private(jsource, jmat, tempmat)
        
        do itarget = 1,nch
          imat = 1 + (itarget-1)*k

          call chunk_neighbors(itarget, wgeo, ibefore, iafter)

          !
          ! for all target chunks, loop over source chunks
          !
          
          do jsource = 1,nch
            jmat = 1 + (jsource-1)*k

            ifoff = 1
            
            if (jsource .eq. ibefore) then
              call znear_buildmat(itarget, jsource, k, wgeo, &
                  fkernel, q1, q2, fgreens, par0, pars1, pars2, &
                  nquad1, xs1, whts1, tempmat)
              call zinsertmat(k, k, tempmat, imat, jmat, &
                  ntot, ntot, zmat)
              ifoff = -1
            endif

            if (jsource .eq. iafter) then
              call znear_buildmat(itarget, jsource, k, wgeo, &
                  fkernel, q1, q2, fgreens, par0, pars1, pars2, &
                  nquad1, xs1, whts1, tempmat)
              call zinsertmat(k, k, tempmat, imat, jmat, &
                  ntot, ntot, zmat)
              ifoff = -1
            endif

            if (jsource .eq. itarget) then
              call zdiagonal_buildmat(itarget, k, wgeo, fkernel, &
                  q1, q2, fgreens, par0, pars1, pars2, xs0, &
                  whts0, tempmat)
              call zinsertmat(k, k, tempmat, imat, jmat, &
                  ntot, ntot, zmat)
              ifoff = -1
            endif

            if (ifoff .eq. 1) then
              call zfar_buildmat(itarget, jsource, k, wgeo, &
                  fkernel, q1, q2, fgreens, par0, pars1, pars2, xs, &
                  whts, tempmat)
              call zinsertmat(k, k, tempmat, imat, jmat, ntot, &
                  ntot, zmat)
            endif

          enddo 
        enddo

        !$omp end parallel do

      end subroutine





      subroutine zfar_buildmat(ich_targ, ich_src, k, wgeo, &
          fkernel, q1, q2, fgreens, par0, pars1, pars2, xs, &
          whts, amat)
        implicit real *8 (a-h,o-z)
        integer, intent(in) :: k
        real *8, intent(in) :: wgeo(*), xs(k), whts(k)
        complex *16, intent(in) :: par0, pars1(*), pars2(*)
        complex *16, intent(in) :: q1(*), q2(*)
        complex *16, intent(out) :: amat(k,k)
        external fkernel, fgreens

        real *8 :: chunk_src(2000), chunk_targ(2000)
        real *8 :: der_src(2000), der_targ(2000)

        !
        ! build a kth order smooth quadrature on for the far
        ! field interation between targ and src chunks
        !

        call chunkaddr(ich_src, wgeo, k_src, ichunk_src, &
             iadj_src, ider_src, ider2_src, h_src)

        call chunkaddr(ich_targ, wgeo, k_targ, ichunk_targ, &
             iadj_targ, ider_targ, ider2_targ, h_targ)

        len = 2*k
        call chunks_dcopy(len, wgeo(ichunk_src), chunk_src)
        call chunks_dcopy(len, wgeo(ider_src), der_src)

        call chunks_dcopy(len, wgeo(ichunk_targ), chunk_targ)
        call chunks_dcopy(len, wgeo(ider_targ), der_targ)

        call zfar_buildmat0(k, chunk_targ, der_targ, &
            h_targ, chunk_src, der_src, h_src, fkernel, q1, q2, &
            fgreens, par0, pars1, pars2, xs, whts, amat)

        return
      end subroutine


      subroutine zfar_buildmat0(k, chunk_targ, &
          der_targ, h_targ, chunk_src, der_src, h_src, fkernel, &
          q1, q2, fgreens, par0, pars1, pars2, xs, whts, amat)
        implicit real *8 (a-h,o-z)
        integer, intent(in) :: k
        real *8, intent(in) :: chunk_targ(2,k), der_targ(2,k), h_targ
        real *8, intent(in) :: chunk_src(2,k), der_src(2,k), h_src
        complex *16, intent(in) :: par0, pars1(*), pars2(*)
        complex *16, intent(in) :: q1(*), q2(*)
        real *8, intent(in) :: xs(k), whts(k)
        complex *16, intent(out) :: amat(k,k)
        external fgreens

        real *8 :: src(10), targ(10)
        real *8 :: target_normal(10), source_normal(10)
        complex *16 :: val

        !
        ! fill up the matrix with the smooth quadrature
        !

        done = 1

        do i = 1,k

          targ(1) = chunk_targ(1,i)
          targ(2) = chunk_targ(2,i)

          xp_targ = der_targ(1,i)
          yp_targ = der_targ(2,i)
          ds_targ = sqrt(xp_targ**2 + yp_targ**2)

          target_normal(1) = yp_targ/ds_targ
          target_normal(2) = -xp_targ/ds_targ

          do j = 1,k

            src(1) = chunk_src(1,j)
            src(2) = chunk_src(2,j)

            xp = der_src(1,j)
            yp = der_src(2,j)

            dd = sqrt(xp**2 + yp**2)
            source_normal(1) = yp/dd
            source_normal(2) = -xp/dd

            call fkernel(par0, src, targ, source_normal, &
                target_normal, q1, q2, fgreens, pars1, pars2, val)

            rnorm = sqrt(xp**2+yp**2)
            dsdt = rnorm*h_src
            amat(i,j) = val*whts(j)*dsdt
          enddo
        enddo

        return
      end subroutine





      subroutine znear_buildmat(ich_targ, ich_src, k, wgeo, &
           fkernel, q1, q2, fgreens, par0, pars1, pars2, nquad, &
           xs, whts, amat)
        implicit real *8 (a-h,o-z)
        integer, intent(in) :: k, nquad
        real *8, intent(in) :: wgeo(*), xs(nquad), whts(nquad)
        complex *16, intent(in) :: par0, pars1(*), pars2(*)
        complex *16, intent(in) :: q1(*), q2(*)
        complex *16, intent(out) :: amat(k,k)
        external fkernel, fgreens

        !
        ! build the matrix performing a neighbor interaction using
        ! bremer style quadratures
        !

        real *8 :: chunk_src(1000)
        real *8 :: chunk_targ(1000), der_src(1000), der_targ(1000)

        call chunkaddr(ich_src, wgeo, k_src, ichunk_src, &
             iadj_src, ider_src, ider2_src, h_src)

        call chunkaddr(ich_targ, wgeo, k_targ, ichunk_targ, &
             iadj_targ, ider_targ, ider2_targ, h_targ)

        len = 2*k
        !!call chunks_dcopy(len, wgeo(ichunk_src), chunk_src)
        !!call chunks_dcopy(len, wgeo(ider_src), der_src)

        !!call chunks_dcopy(len, wgeo(ichunk_targ), chunk_targ)
        !!call chunks_dcopy(len, wgeo(ider_targ), der_targ)

        call znear_buildmat0(k, wgeo(ichunk_targ), wgeo(ider_targ), &
            h_targ, wgeo(ichunk_src), wgeo(ider_src), h_src, &
            fkernel, q1, q2, fgreens, par0, pars1, pars2, nquad, &
            xs, whts, amat)

        !!call znear_buildmat0(k, wgeo(ichunk_targ)chunk_targ, der_targ, &
        !!    h_targ, chunk_src, der_src, h_src, fkernel, q1, q2, &
        !!    fgreens, par0, pars1, pars2, nquad, xs, whts, amat)

        return
      end subroutine


      subroutine znear_buildmat0(k, chunk_targ, der_targ, &
            h_targ, chunk_src, der_src, h_src, fkernel, q1, q2, &
            fgreens, par0, pars1, pars2, nquad, xs, whts, amat)
        implicit real *8 (a-h,o-z)
        integer, intent(in) :: k, nquad
        real *8, intent(in) :: chunk_targ(2,k), der_targ(2,k), h_targ
        real *8, intent(in) :: chunk_src(2,k), der_src(2,k), h_src
        complex *16, intent(in) :: par0, pars1(*), pars2(*)
        complex *16, intent(in) :: q1(*), q2(*)
        real *8, intent(in) :: xs(nquad), whts(nquad)
        complex *16, intent(out) :: amat(k,k)
        external fgreens

        real *8 :: src(10), targ(10), source_normal(10)
        real *8 :: ts(1000), ainterp(1000000), target_normal(10)
        real *8 :: xcoefs_src(1000), ycoefs_src(1000)
        real *8 :: xpcoefs_src(1000), ypcoefs_src(1000)
        complex *16 :: val, temp(10000), w(1000000)

        !
        ! build the k by k near interaction matrix.
        !

        done = 1
        call chunksexps(k, chunk_src, xcoefs_src, ycoefs_src)
        call chunksexps(k, der_src, xpcoefs_src, ypcoefs_src)

        call lematrin(k, nquad, xs, ainterp, ts, w)

        do i = 1,k

          targ(1) = chunk_targ(1,i)
          targ(2) = chunk_targ(2,i)

          xp_targ = der_targ(1,i)
          yp_targ = der_targ(2,i)
          ds_targ = sqrt(xp_targ**2 + yp_targ**2)

          target_normal(1) = yp_targ/ds_targ
          target_normal(2) = -xp_targ/ds_targ


          do j = 1,nquad

            t_src = xs(j)
            call legeexev(t_src, src(1), xcoefs_src, k-1)
            call legeexev(t_src, src(2), ycoefs_src, k-1)
            call legeexev(t_src, xp, xpcoefs_src, k-1)
            call legeexev(t_src, yp, ypcoefs_src, k-1)

            dd = sqrt(xp**2 + yp**2)
            source_normal(1) = yp/dd
            source_normal(2) = -xp/dd

            call fkernel(par0, src, targ, source_normal, &
                target_normal, q1, q2, fgreens, pars1, pars2, val)

            dsdt = sqrt(xp**2+yp**2)*h_src
            temp(j) = val*whts(j)*dsdt
          enddo

          !!!!call lematrin(k, nquad, xs, ainterp, ts, w)
          i1 = 1
          call zmat_by_dmat(i1, nquad, temp, k, ainterp, w)

          do j = 1,k
            amat(i,j) = w(j)
          enddo
        enddo

        return
      end subroutine





      subroutine znear_apply_fast(ich_targ, ich_src, k, wgeo, &
           fkernel, q1, q2, fgreens, par0, pars1, pars2, nquad, &
           xs, whts, ainterp, charge, pot)
        implicit real *8 (a-h,o-z)
        integer, intent(in) :: k, nquad
        real *8, intent(in) :: wgeo(*), xs(nquad), whts(nquad)
        real *8, intent(in) :: ainterp(nquad,k)
        complex *16, intent(in) :: par0, pars1(*), pars2(*)
        complex *16, intent(in) :: q1(*), q2(*), charge(k)
        complex *16, intent(out) :: pot(k)
        external fkernel, fgreens

        real *8 :: xs7(10000), ts(100000), w(1000000)
        !
        ! build the matrix performing a neighbor interaction using
        ! bremer style quadratures
        !

        call chunkaddr(ich_src, wgeo, k_src, ichunk_src, &
             iadj_src, ider_src, ider2_src, h_src)

        call chunkaddr(ich_targ, wgeo, k_targ, ichunk_targ, &
             iadj_targ, ider_targ, ider2_targ, h_targ)

        !!!!call lematrin(k, nquad, xs, ainterp, ts, w)

        !!!!call prin2('ainterp = *', ainterp, 30)

        call znear_apply_fast0(k, wgeo(ichunk_targ), wgeo(ider_targ), &
            h_targ, wgeo(ichunk_src), wgeo(ider_src), h_src, &
            fkernel, q1, q2, fgreens, par0, pars1, pars2, nquad, &
            xs, whts, ainterp, charge, pot)

        return
      end subroutine


      subroutine znear_apply_fast0(k, chunk_targ, der_targ, &
            h_targ, chunk_src, der_src, h_src, fkernel, q1, q2, &
            fgreens, par0, pars1, pars2, nquad, xs, whts, ainterp, &
            charge, pot)
        implicit real *8 (a-h,o-z)
        integer, intent(in) :: k, nquad, ainterp(nquad,k)
        real *8, intent(in) :: chunk_targ(2,k), der_targ(2,k), h_targ
        real *8, intent(in) :: chunk_src(2,k), der_src(2,k), h_src
        complex *16, intent(in) :: par0, pars1(*), pars2(*)
        complex *16, intent(in) :: q1(*), q2(*), charge(k)
        real *8, intent(in) :: xs(nquad), whts(nquad)
        complex *16, intent(out) :: pot(k)
        external fgreens

        real *8 :: src(10), targ(10), source_normal(10)
        real *8 :: target_normal(10)
        real *8 :: xcoefs_src(1000), ycoefs_src(1000)
        real *8 :: xpcoefs_src(1000), ypcoefs_src(1000)
        real *8 :: xys(2,10000), dxys(2,10000)
        complex *16 :: val, ts(100000), w(1000000)
        complex *16 :: chargenew(10000)

        !
        ! merely apply the k by k near interaction matrix.
        !

        done = 1
        !!!!call lematrin(k, nquad, xs, ainterp, ts, w)

        
        i1 = 1
        call dmat_by_zmat(nquad, k, ainterp, i1, charge, chargenew)
        call mat2vec(nquad, k, ainterp, chunk_src, xys)
        call mat2vec(nquad, k, ainterp, der_src, dxys)
        
        do i = 1,k

          targ(1) = chunk_targ(1,i)
          targ(2) = chunk_targ(2,i)

          xp_targ = der_targ(1,i)
          yp_targ = der_targ(2,i)
          ds_targ = sqrt(xp_targ**2 + yp_targ**2)

          target_normal(1) = yp_targ/ds_targ
          target_normal(2) = -xp_targ/ds_targ

          pot(i) = 0
          
          do j = 1,nquad

            src(1) = xys(1,j)
            src(2) = xys(2,j)
            xp = dxys(1,j)
            yp = dxys(2,j)
            
            dd = sqrt(xp**2 + yp**2)
            source_normal(1) = yp/dd
            source_normal(2) = -xp/dd

            call fkernel(par0, src, targ, source_normal, &
                target_normal, q1, q2, fgreens, pars1, pars2, val)

            dsdt = sqrt(xp**2+yp**2)*h_src
            pot(i) = pot(i) + val*whts(j)*dsdt*chargenew(j)
          enddo

        enddo

        return
      end subroutine





      subroutine mat2vec(m, n, a, xs, ys)
        implicit real *8 (a-h,o-z)
        real *8 :: a(m,n), xs(2,n), ys(2,m)

        do i = 1,m
          dd1 = 0
          dd2 = 0
          do j = 1,n
            dd1 = dd1 + a(i,j)*xs(1,j)
            dd2 = dd2 + a(i,j)*xs(2,j)
          enddo
          ys(1,i) = dd1
          ys(2,i) = dd2
        enddo

        return
      end subroutine mat2vec




      
      subroutine zdiagonal_buildmat(ich, k, wgeo, fkernel, q1, q2, &
          fgreens, par0, pars1, pars2, xs0, whts0, amat)
        implicit real *8 (a-h,o-z)
        integer, intent(in) :: k
        real *8, intent(in) :: wgeo(*), xs0(k,k), whts0(k,k)
        complex *16, intent(in) :: par0, pars1(*), pars2(*)
        complex *16, intent(in) :: q1(*), q2(*)
        complex *16, intent(out) :: amat(k,k)
        external fkernel, fgreens

        !
        ! build the k by k self interaction matrix. this routine (will) builds
        ! matrices which are square-root weights by default, since it's
        ! the morally correct thing to do. this uses Bremer quadratures.
        !
        ! input:
        !   ich - the chunk number
        !   k - dimension of amat, i.e. the number of points on the chunk
        !   wgeo - geometry structure created by chunks.f
        !   fkernel - the arbitrary kernel routine, see above
        !       documention since this routine receives the greens function
        !   fgreens - the single layer greens function routine
        !       it must have the calling sequence
        !           fgreens(par0, src, targ, pars1, pars2, val)
        !       and val is assumed to be real
        !   par0, pars1, pars2 - inputs to fgreens
        !
        ! output:
        !   amat - the self-interaction matrix
        !

        call chunkaddr(ich, wgeo, k, ichunk, iadj, ider, ider2, h)

        call zdiagonal_buildmat0(k, wgeo(ichunk), &
            wgeo(ider), h, fkernel, q1, q2, fgreens, par0, &
            pars1, pars2, xs0, whts0, amat)

        return
      end subroutine


      subroutine zdiagonal_buildmat0(k, xys, dxys, h, &
          fkernel, q1, q2, fgreens, par0, pars1, pars2, &
          xs, whts, amat)
        implicit real *8 (a-h,o-z)
        integer, intent(in) :: k 
        real *8, intent(in) :: xys(2,k), dxys(2,k), h
        complex *16, intent(in) :: par0, pars1(*), pars2(*)
        complex *16, intent(in) :: q1(*), q2(*)
        real *8, intent(in) :: xs(k,k), whts(k,k)
        complex *16, intent(out) :: amat(k,k)
        external fgreens

        real *8 :: src(10), targ(10), source_normal(10)
        real *8 :: ts(1000), ainterp(1000000), target_normal(10)
        real *8 :: xcoefs(1000)
        real *8 :: ycoefs(1000), xpcoefs(1000), ypcoefs(1000)
        complex *16 :: val, temp(10000), w(1000000)

        !
        ! build the k by k self interaction matrix. no square-root
        ! weighting yet... that is handled in the outer routine,
        ! xxx_xxx_direct.
        ! Note: k should not be larger than 800
        !

        done = 1
        call chunksexps(k, xys, xcoefs, ycoefs)
        call chunksexps(k, dxys, xpcoefs, ypcoefs)

        do i = 1,k

          targ(1) = xys(1,i)
          targ(2) = xys(2,i)

          xp_targ = dxys(1,i)
          yp_targ = dxys(2,i)
          ds_targ = sqrt(xp_targ**2 + yp_targ**2)

          target_normal(1) = yp_targ/ds_targ
          target_normal(2) = -xp_targ/ds_targ

          do j = 1,k

            t_src = xs(j,i)
            call legeexev(t_src, src(1), xcoefs, k-1)
            call legeexev(t_src, src(2), ycoefs, k-1)
            call legeexev(t_src, xp, xpcoefs, k-1)
            call legeexev(t_src, yp, ypcoefs, k-1)

            dd = sqrt(xp**2 + yp**2)
            source_normal(1) = yp/dd
            source_normal(2) = -xp/dd

            call fkernel(par0, src, targ, source_normal, &
                target_normal, q1, q2, fgreens, pars1, pars2, val)

            dsdt = sqrt(xp**2+yp**2)*h
            temp(j) = val*whts(j,i)*dsdt
          enddo

          call lematrin(k, k, xs(1,i), ainterp, ts, w)
          i1 = 1
          call zmat_by_dmat(i1, k, temp, k, ainterp, w)

          do j = 1,k
            amat(i,j) = w(j)
          enddo
        enddo

        return
      end subroutine





      subroutine zdiagonal_apply_fast(k, xys, dxys, h, &
          fkernel, q1, q2, fgreens, par0, pars1, pars2, &
          xs, whts, amat)
        implicit real *8 (a-h,o-z)
        integer, intent(in) :: k 
        real *8, intent(in) :: xys(2,k), dxys(2,k), h
        complex *16, intent(in) :: par0, pars1(*), pars2(*)
        complex *16, intent(in) :: q1(*), q2(*)
        real *8, intent(in) :: xs(k,k), whts(k,k)
        complex *16, intent(out) :: amat(k,k)
        external fgreens

        real *8 :: src(10), targ(10), source_normal(10)
        real *8 :: ts(1000), ainterp(1000000), target_normal(10)
        real *8 :: xcoefs(1000)
        real *8 :: ycoefs(1000), xpcoefs(1000), ypcoefs(1000)
        complex *16 :: val, temp(10000), w(1000000)

        !
        ! build the k by k self interaction matrix. no square-root
        ! weighting yet... that is handled in the outer routine,
        ! xxx_xxx_direct.
        ! Note: k should not be larger than 800
        !

        done = 1
        call chunksexps(k, xys, xcoefs, ycoefs)
        call chunksexps(k, dxys, xpcoefs, ypcoefs)

        do i = 1,k

          targ(1) = xys(1,i)
          targ(2) = xys(2,i)

          xp_targ = dxys(1,i)
          yp_targ = dxys(2,i)
          ds_targ = sqrt(xp_targ**2 + yp_targ**2)

          target_normal(1) = yp_targ/ds_targ
          target_normal(2) = -xp_targ/ds_targ

          
          do j = 1,k

            t_src = xs(j,i)
            call legeexev(t_src, src(1), xcoefs, k-1)
            call legeexev(t_src, src(2), ycoefs, k-1)
            call legeexev(t_src, xp, xpcoefs, k-1)
            call legeexev(t_src, yp, ypcoefs, k-1)

            dd = sqrt(xp**2 + yp**2)
            source_normal(1) = yp/dd
            source_normal(2) = -xp/dd

            call fkernel(par0, src, targ, source_normal, &
                target_normal, q1, q2, fgreens, pars1, pars2, val)

            dsdt = sqrt(xp**2+yp**2)*h
            temp(j) = val*whts(j,i)*dsdt
          enddo

          call lematrin(k, k, xs(1,i), ainterp, ts, w)

          call dmat_by_zvec(k, k, ainterp, charge, w)

          !!i1 = 1
          !!call zmat_by_dmat(i1, k, temp, k, ainterp, w)

          !!!!!!call zvecvec(k, temp, w, pot(i))
          
          !!do j = 1,k
          !!  amat(i,j) = w(j)
          !!enddo
        enddo

        return
      end subroutine zdiagonal_apply_fast





      subroutine zvecvec(n, x, y, prod)
        implicit real *8 (a-h,o-z)
        complex *16 :: x(n), y(n), prod

        prod = 0
        do i = 1,n
          prod = prod + x(i)*y(i)
        enddo

        return
      end subroutine zvecvec




      
      subroutine dmat_by_zvec(m, n, a, x, y)
        implicit real *8 (a-h,o-z)
        real *8 :: a(m,n)
        complex *16 :: x(m), y(n), cd

        do i = 1,m
          cd = 0
          do j = 1,n
            cd = cd + a(i,j)*x(j)
          enddo
          y(i) = cd
        enddo

        return
      end subroutine dmat_by_zvec




      
      subroutine zinsertmat(km, kn, amat, iloc, jloc, m, n, cmat)
        implicit real *8 (a-h,o-z)
        complex *16 :: amat(km,kn), cmat(m,n)

        do j = 1,kn
          do i = 1,km
            cmat(iloc-1+i,jloc-1+j) = amat(i,j)
          enddo
        enddo

        return
      end subroutine

