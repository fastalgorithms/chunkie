
      subroutine multichunkfuncrmodes(eps, ifclosed, chsmall, tas, tbs,
     1     pars, iparsindeces, k, novers, wgeos, lwgeos,
     2     ncomp, ncompmax, lused, ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     this subroutine forms boundary components specified by funcurve
c     and various parameters and stores their discretization in 
c     the array wgeos as a multichunk.
c
c     INPUT
c
c     eps - real *8 array. eps(i) is the tolerance passed to chunkfunc
c           for the ith boundary component
c     ifclosed - integer array. ifclosed(i) = 1 if the ith component is
c           to be treated as a closed curve by chunkfunc, 
c           ifclosed(i) = 0 otherwise.
c     chsmall - real *8 array. chsmall(i) is the size of the smallest
c               chunk near the endpoints of the ith component, if the
c               ith component is to be treated as a closed curve.
c     tas, tbs - real *8 array. the start and end points of the 
c                ith component in parameter space are tas(i), tbs(i)
c
c     TO DO: add docs for these pars and iparsindeces... designed to
c     be called by funcurve_by_mode_cmextra(t,w,x,y,dxdt,dydt,dxdt2,dydt2)
c
c     pars - real *8 array of parameters.
c     iparsindeces - integer array. iparsindeces(i) determines
c                    where in the large pars array the parameters
c                    for the ith component begin. That is, when
c                    the ith chunk is discretized, the subroutine
c                    funcurve is called with
c
c         funcurve(t,pars(iparsindeces(i)),x,y,dxdt,dydt,dxdt2,dydt2)
c
c     k - integer. the order of discretization to use on the chunks
c     novers - integer array. novers(i) is the oversampling factor
c              to be used on the ith chunk
c     lwgeos - integer. the length of the wgeos array
c     ncomp - the number of boundary components (the eps, chsmall,
c             etc. arrays should be length ncomp)
c     ncompmax - the maximum number of components to be stored in the
c                multichunk (wgeos is initialized with this setting)
c                note, it is necessary that ncompmax >= ncomp
c
c     OUTPUT
c
c     wgeos - the multichunk for these boundary components
c     lused - the length of the wgeos array which was used.
c     ier - integer, flag. error flag.
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real *8 eps(*), chsmall(*), tas(*), tbs(*), pars(*)
      real *8 wgeos(*)
      integer ncomp, ncompmax, iparsindeces(*), k, novers(*)
      integer lused, ier, lwgeos, ifclosed(*)
c     local variables
      external funcurve_by_mode_cmextra
      

      call multichunkfunc(eps, ifclosed, chsmall, tas, tbs,
     1     funcurve_by_mode_cmextra, pars, iparsindeces, k, novers,
     2     wgeos, lwgeos, ncomp, ncompmax, lused, ier)
      
      return
      end
      


      subroutine funcurve_by_mode_cmextra(t,w,x,y,dxdt,dydt,dxdt2,dydt2)
      implicit real *8 (a-h,o-z)
      real *8 t, w(*), x, y, dxdt, dydt, dxdt2, dydt2
      real *8 cx, cy, r, rp, rpp, dxdr, dydr, dxdtheta, dydtheta
      real *8 tt, pi2

      pi2 = 8.0d0*datan(1.0d0)

      nmodes = w(1)
      imodes = w(2)
      iort = w(3)
      cx = w(4)
      cy = w(5)

      tt = t
      if (iort .eq. -1) then
         tt = pi2-t
      endif


      index = imodes

      r = w(index)
      index = index + 1
      rp = 0.0d0
      rpp = 0.0d0

      do i = 1,nmodes
         temp = w(index)
         index = index+1
         r = r+temp*dcos(i*tt)
         rp = rp-i*temp*dsin(i*tt)
         rpp = rpp-i**2*temp*dcos(i*tt)
         temp = w(index)
         index = index+1
         r = r+temp*dsin(i*tt)
         rp = rp+i*temp*dcos(i*tt)
         rpp = rpp-i**2*temp*dsin(i*tt)
      enddo

      x = cx+r*dcos(tt)
      y = cy+r*dsin(tt)

      dxdr = dcos(tt)
      dydr = dsin(tt)
        
      dxdtheta = -r*dsin(tt)
      dxdtheta2 = -r*cos(tt)
      dydtheta = r*dcos(tt)
      dydtheta2 = -r*dsin(tt)

      dxdrdtheta = -dsin(tt)
      dydrdtheta = dcos(tt)

      dxdt = iort*(dxdr*rp+dxdtheta)
      dydt = iort*(dydr*rp+dydtheta)

      dxdt2 = dxdr*rpp+dxdrdtheta*2.0d0*rp+dxdtheta2
      dydt2 = dydr*rpp+dydrdtheta*2.0d0*rp+dydtheta2

      return
      end

      subroutine multichunkfunc_file_by_mode(ifile,wgeos,
     1     lwgeos,ncompmax,ncomp,k,lused,ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      integer ifile, lwgeos, ncompmax, ier
      real *8 wgeos(*)
      integer ncomp, k, lused
c     local variables
      real *8, allocatable :: eps(:), chsmall(:), tas(:), tbs(:)
      real *8, allocatable :: pars(:)
      integer, allocatable :: ifclosed(:), iparsindeces(:), novers(:)
      integer lpars
      external funcurve_by_mode_cmextra

      ier = 0

c     initialize multichunk

      call multichunk_init(wgeos,lwgeos,ncompmax,ier)

c     find out number of components and length of pars array

      call query_multichunk_file_by_mode(ifile,
     1     ncomp, lpars, ier)

      allocate(eps(ncomp))
      allocate(chsmall(ncomp))
      allocate(tas(ncomp))
      allocate(tbs(ncomp))
      allocate(pars(lpars))
      allocate(ifclosed(ncomp))
      allocate(iparsindeces(ncomp))
      allocate(novers(ncomp))

c     read in curve data (used by funcurve_by_mode_cmextra)

      call read_multichunk_file_by_mode(ifile,
     1     ncomp, eps, ifclosed, chsmall, tas, tbs, pars,
     2     iparsindeces, k, novers, ier)

c     make multichunk for this geometry

      call multichunkfunc(eps, ifclosed, chsmall, tas, tbs,
     1     funcurve_by_mode_cmextra, pars, iparsindeces, k, novers, 
     2     wgeos, lwgeos, ncomp, ncompmax, lused, ier)

      return
      end
      

      subroutine read_multichunk_file_by_mode(ifile,
     1     ncomp, eps, ifclosed, chsmall, tas, tbs, pars,
     2     iparsindeces, k, novers, ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INTEGER
c
c     INPUT 
c
c     ifile - integer for file stream
c     
c     OUTPUT
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      integer ifile, ncomp, ifclosed(*), iparsindeces(*), k
      integer novers(*), ier
      real *8 eps(*), chsmall(*), tas(*), tbs(*), pars(*)
c     local variables
      integer i, icomp, iparstemp, j, nmodes, ios, iort
      real *8 coefftemp, pi, cx, cy, coscoefftemp, sincoefftemp

      pi = 4.0d0*datan(1.0d0)
      
      ier = 0
      iparstemp = 1
      rewind(ifile,iostat=ios)
      if (ios .ne. 0) write(*,*) 'BAD REWIND'
      read(ifile,*,iostat=ios) ncomp
      if (ios .ne. 0) goto 1234
      read(ifile,*,iostat=ios) k
      if (ios .ne. 0) goto 1234
     
      do i = 1,ncomp

c     first should be component number

         read(ifile,*,iostat=ios) icomp
         if (ios .ne. 0) goto 1234
         if (i .ne. icomp) then
            write(*,*) 'ERROR READING FILE, WRONG FORMAT: '
            write(*,*) 'Wrong component number'
            ier = 1
            return
         endif

c     because it's described by modes, some
c     parameters are automatic

         ifclosed(i) = 1
         chsmall(i) = 1.0d0
         tas(i) = 0.0d0
         tbs(i) = 2.0d0*pi

c     get info for chunkfunc routine and orientation

         read(ifile,*,iostat=ios) eps(i)
         if (ios .ne. 0) goto 1234
         read(ifile,*,iostat=ios) iort
         if (ios .ne. 0) goto 1234
         read(ifile,*,iostat=ios) novers(i)
         if (ios .ne. 0) goto 1234

c     then parameters passed to chunkfunc routine (info
c     about curve parameterization)

         iparsindeces(i) = iparstemp
         
c     first is number of modes

         read(ifile,*,iostat=ios) nmodes
         if (ios .ne. 0) goto 1234

         if (nmodes .lt. 0) then
            ier = 2
            write(*,*) 'ERROR READING FILE: '
            write(*,*) 'NMODES negative for component', i
            return
         endif

         pars(iparstemp) = nmodes+0.1d0
         iparstemp = iparstemp+1

         pars(iparstemp) = 6+0.1d0
         iparstemp = iparstemp+1

         pars(iparstemp) = iort*(1.1d0)
         iparstemp = iparstemp+1

c     then center of parameterized curve

         read(ifile,*,iostat=ios) cx, cy
         if (ios .ne. 0) goto 1234
         pars(iparstemp) = cx
         iparstemp = iparstemp+1
         pars(iparstemp) = cy
         iparstemp = iparstemp+1

c     then the modes (nmodes should be >= 0
c     and there should be a zero-th mode)
c     because it's a cosine and sine series,
c     there should be 2*nmodes+1 coefficients
c     describing the curve

c     zero-th mode (constant term). MUST EXIST
         
         read(ifile,*,iostat=ios) coefftemp
         if (ios .ne. 0) goto 1234
         pars(iparstemp) = coefftemp
         iparstemp = iparstemp+1
         
c     higher modes

         do j = 1,nmodes
            
c     cos(jx), sin(jx) coefficients

            read(ifile,*,iostat=ios) coscoefftemp, sincoefftemp
            if (ios .ne. 0) goto 1234
            pars(iparstemp) = coscoefftemp
            iparstemp = iparstemp+1
            pars(iparstemp) = sincoefftemp
            iparstemp = iparstemp+1         
         enddo

      enddo
         

      return

 1234 continue

      ier = 3
      write(*,*) 'ERROR READING FILE:'
      write(*,*) 'REACHED END OR BAD READ'
      write(*,*) 'CHECK FILE FORMAT'
      return

      end


      subroutine query_multichunk_file_by_mode(ifile,
     1     ncomp, lpars, ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INTEGER
c
c     INPUT 
c
c     ifile - integer for file stream
c     
c     OUTPUT
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      integer ifile, ncomp
      integer ier, lpars
c     local variables
      integer i, icomp, iparstemp, j, nmodes, ios, k
      real *8 coefftemp, pi, cx, cy, coscoefftemp, sincoefftemp
      real *8 temp

      pi = 4.0d0*datan(1.0d0)
      
      ier = 0
      iparstemp = 1
      rewind(ifile,iostat=ios)
      if (ios .ne. 0) write(*,*) 'BAD REWIND'
      read(ifile,*,iostat=ios) ncomp
      if (ios .ne. 0) goto 1235
      read(ifile,*,iostat=ios) k
      if (ios .ne. 0) goto 1235
     
      do i = 1,ncomp

c     first should be component number

         read(ifile,*,iostat=ios) icomp
         if (ios .ne. 0) goto 1235
         if (i .ne. icomp) then
            write(*,*) 'ERROR READING FILE, WRONG FORMAT: '
            write(*,*) 'Wrong component number'
            ier = 1
            return
         endif

c     get info for chunkfunc routine 

         read(ifile,*,iostat=ios) temp
         if (ios .ne. 0) goto 1235
         read(ifile,*,iostat=ios) temp
         if (ios .ne. 0) goto 1235
         read(ifile,*,iostat=ios) temp
         if (ios .ne. 0) goto 1235

c     then parameters passed to chunkfunc routine (info
c     about curve parameterization)

c     first is number of modes

         read(ifile,*,iostat=ios) nmodes
         if (ios .ne. 0) goto 1235

         if (nmodes .lt. 0) then
            ier = 2
            write(*,*) 'ERROR READING FILE: '
            write(*,*) 'NMODES negative for component', i
            return
         endif

         iparstemp = iparstemp+1
         iparstemp = iparstemp+1
         iparstemp = iparstemp+1



c     then center of parameterized curve

         read(ifile,*,iostat=ios) cx, cy
         if (ios .ne. 0) goto 1235
         iparstemp = iparstemp+1
         iparstemp = iparstemp+1

c     then the modes (nmodes should be >= 0
c     and there should be a zero-th mode)
c     because it's a cosine and sine series,
c     there should be 2*nmodes+1 coefficients
c     describing the curve

c     zero-th mode (constant term). MUST EXIST
         
         read(ifile,*,iostat=ios) coefftemp
         if (ios .ne. 0) goto 1235
         iparstemp = iparstemp+1
         
c     higher modes

         do j = 1,nmodes
            
c     cos(jx), sin(jx) coefficients

            read(ifile,*,iostat=ios) coscoefftemp, sincoefftemp
            if (ios .ne. 0) goto 1235
            iparstemp = iparstemp+1
            iparstemp = iparstemp+1         
         enddo

      enddo

      lpars = iparstemp-1

      return

 1235 continue

      ier = 3
      write(*,*) 'ERROR READING FILE:'
      write(*,*) 'REACHED END OR BAD READ'
      write(*,*) 'CHECK FILE FORMAT'
      return

      end


      subroutine funcurve_fourierxy_cmextra(t,w,x,y,dxdt,dydt,dxdt2,
     1     dydt2)
      implicit real *8 (a-h,o-z)
      real *8 t, w(*), x, y, dxdt, dydt, dxdt2, dydt2
      real *8 cx, cy, r, rp, rpp, dxdr, dydr, dxdtheta, dydtheta
      real *8 tt, pi2
      complex *16 xcoeff, ycoeff, wi, eye, zmul
      data eye /(0.0d0,1.0d0) /
      pi2 = 8.0d0*datan(1.0d0)

      nmodes = w(1)
      imodesx = w(2)
      imodesy = w(3)

      xcoeff = dcmplx(w(imodesx),w(imodesx+1))
      ycoeff = dcmplx(w(imodesx),w(imodesx+1))
      
      x = dreal(xcoeff)
      y = dreal(ycoeff)
      dxdt = 0.0d0
      dydt = 0.0d0
      dxdt2 = 0.0d0
      dydt2 = 0.0d0

      wi = exp(eye*pi2*t)
      zmul = wi

      do i = 1,nmodes
         xcoeff = dcmplx(w(imodesx+2*i),w(imodesx+2*i+1))
         ycoeff = dcmplx(w(imodesy+2*i),w(imodesy+2*i+1))
         x = x+2.0d0*dreal(xcoeff*wi)
         y = y+2.0d0*dreal(ycoeff*wi)
         dxdt = dxdt + 2.0d0*pi2*i*dreal(xcoeff*wi*eye)
         dydt = dydt + 2.0d0*pi2*i*dreal(ycoeff*wi*eye)
         dxdt2 = dxdt2 - 2.0d0*pi2**2*i**2*dreal(xcoeff*wi)
         dydt2 = dydt2 - 2.0d0*pi2**2*i**2*dreal(ycoeff*wi)
         wi = wi*zmul
      enddo

      return
      end

      subroutine multichunkfunc_file_fourierxy(ifile,wgeos,
     1     lwgeos,ncompmax,ncomp,k,lused,ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      integer ifile, lwgeos, ncompmax, ier
      real *8 wgeos(*)
      integer ncomp, k, lused
c     local variables
      real *8, allocatable :: eps(:), chsmall(:), tas(:), tbs(:)
      real *8, allocatable :: pars(:)
      integer, allocatable :: ifclosed(:), iparsindeces(:), novers(:)
      integer, allocatable :: iort(:)
      integer lpars, i
      external funcurve_fourierxy_cmextra

c     initialize multichunk

      call multichunk_init(wgeos,lwgeos,ncompmax,ier)

c     find out number of components and length of pars array

      call query_multichunk_file_fourierxy(ifile,
     1     ncomp, lpars, ier)

      allocate(eps(ncomp))
      allocate(chsmall(ncomp))
      allocate(tas(ncomp))
      allocate(tbs(ncomp))
      allocate(pars(lpars))
      allocate(ifclosed(ncomp))
      allocate(iparsindeces(ncomp))
      allocate(novers(ncomp))
      allocate(iort(ncomp))

c     read in curve data (used by funcurve_by_mode_cmextra)

      call read_multichunk_file_fourierxy(ifile,
     1     ncomp, eps, ifclosed, chsmall, tas, tbs, pars,
     2     iparsindeces, iort, k, novers, ier)

c      call prin2('pars *',pars,lpars)
c      call prinf('iparsindeces *',iparsindeces,ncomp)
c      call prin2('eps *',eps,ncomp)
c      call prinf('ifclosed *',ifclosed,ncomp)
c      call prin2('chsmall *',chsmall,ncomp)
c      call prin2('tas *',tas,ncomp)
c      call prin2('tbs *',tbs,ncomp)
c      call prinf('iort *',iort,ncomp)
c      call prinf('k *',k,1)

c     make multichunk for this geometry

      call multichunkfunc(eps, ifclosed, chsmall, tas, tbs,
     1     funcurve_fourierxy_cmextra, pars, iparsindeces, k, novers, 
     2     wgeos, lwgeos, ncomp, ncompmax, lused, ier)

c     reverse components if requested (if iort(icomp) = -1)

      do i = 1,ncomp
         if (iort(i) .eq. -1) then
            call multichunk_reverse_comp(wgeos,i,ier)
            call multichunk_sort_comp(wgeos,i,ier)            
         endif
      enddo


      return
      end
      

      subroutine read_multichunk_file_fourierxy(ifile,
     1     ncomp, eps, ifclosed, chsmall, tas, tbs, pars,
     2     iparsindeces, iort, k, novers, ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INTEGER
c
c     INPUT 
c
c     ifile - integer for file stream
c     
c     OUTPUT
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      integer ifile, ncomp, ifclosed(*), iparsindeces(*), k
      integer novers(*), ier, iort(*)
      real *8 eps(*), chsmall(*), tas(*), tbs(*), pars(*)
c     local variables
      integer i, icomp, iparstemp, j, nmodes, ios, iorttemp
      real *8 coefftemp, pi, cx, cy, realcoefftemp, imagcoefftemp

      pi = 4.0d0*datan(1.0d0)
      
      ier = 0
      iparstemp = 1
      rewind(ifile,iostat=ios)
      if (ios .ne. 0) write(*,*) 'BAD REWIND'
      read(ifile,*,iostat=ios) ncomp
      if (ios .ne. 0) goto 1236
      read(ifile,*,iostat=ios) k
      if (ios .ne. 0) goto 1236
     
      do i = 1,ncomp

c     first should be component number

         read(ifile,*,iostat=ios) icomp
         if (ios .ne. 0) goto 1236
         if (i .ne. icomp) then
            write(*,*) 'ERROR READING FILE, WRONG FORMAT: '
            write(*,*) 'Wrong component number'
            ier = 1
            return
         endif

c     because it's described by modes, some
c     parameters are automatic

         ifclosed(i) = 1
         chsmall(i) = 1.0d0
         tas(i) = 0.0d0
         tbs(i) = 1.0d0

c     get info for chunkfunc routine and orientation

         read(ifile,*,iostat=ios) eps(i)
         if (ios .ne. 0) goto 1236
         read(ifile,*,iostat=ios) iorttemp
         if (ios .ne. 0) goto 1236
         read(ifile,*,iostat=ios) novers(i)
         if (ios .ne. 0) goto 1236

         iort(i) = iorttemp

c     then parameters passed to chunkfunc routine (info
c     about curve parameterization)

         iparsindeces(i) = iparstemp
         
c     first is number of modes

         read(ifile,*,iostat=ios) nmodes
         if (ios .ne. 0) goto 1236

         if (nmodes .lt. 0) then
            ier = 2
            write(*,*) 'ERROR READING FILE: '
            write(*,*) 'NMODES negative for component', i
            return
         endif

         pars(iparstemp) = nmodes+0.1d0
         iparstemp = iparstemp+1

         pars(iparstemp) = 4+0.1d0
         iparstemp = iparstemp+1
         pars(iparstemp) = 4+0.1d0+2*(nmodes+1)
         iparstemp = iparstemp+1

c     then the modes (nmodes should be >= 0
c     and there should be a zero-th mode)
c     coefficients are complex, in columns of
c     two real numbers. First the x coeffs
c     are listed, then the y coeffs

c     x coeffs

         do j = 0,nmodes
            
c     real, imaginary parts of coefficients

            read(ifile,*,iostat=ios) realcoefftemp, imagcoefftemp
            if (ios .ne. 0) goto 1236
            pars(iparstemp) = realcoefftemp
            iparstemp = iparstemp+1
            pars(iparstemp) = imagcoefftemp
            iparstemp = iparstemp+1         
         enddo

c     y coeffs

         do j = 0,nmodes
            
c     real, imaginary parts of coefficients

            read(ifile,*,iostat=ios) realcoefftemp, imagcoefftemp
            if (ios .ne. 0) goto 1236
            pars(iparstemp) = realcoefftemp
            iparstemp = iparstemp+1
            pars(iparstemp) = imagcoefftemp
            iparstemp = iparstemp+1         
         enddo

      enddo
         

      return

 1236 continue

      ier = 3
      write(*,*) 'ERROR READING FILE:'
      write(*,*) 'REACHED END OR BAD READ'
      write(*,*) 'CHECK FILE FORMAT'
      return

      end


      subroutine query_multichunk_file_fourierxy(ifile,
     1     ncomp, lpars, ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INTEGER
c
c     INPUT 
c
c     ifile - integer for file stream
c     
c     OUTPUT
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      integer ifile, ncomp, lpars
      integer ier
c     local variables
      integer i, icomp, itemp, nmodes, ios, iparstemp, j
      real *8 dtemp, dtemp2
      
      ier = 0
      iparstemp = 1
      rewind(ifile,iostat=ios)
      if (ios .ne. 0) write(*,*) 'BAD REWIND'
      read(ifile,*,iostat=ios) ncomp
      if (ios .ne. 0) goto 1237
      read(ifile,*,iostat=ios) itemp
      if (ios .ne. 0) goto 1237
     
      do i = 1,ncomp

c     first should be component number

         read(ifile,*,iostat=ios) icomp
         if (ios .ne. 0) goto 1237
         if (i .ne. icomp) then
            write(*,*) 'ERROR READING FILE, WRONG FORMAT: '
            write(*,*) 'Wrong component number'
            ier = 1
            return
         endif

c     get info for chunkfunc routine and orientation

         read(ifile,*,iostat=ios) dtemp
         if (ios .ne. 0) goto 1237
         read(ifile,*,iostat=ios) itemp
         if (ios .ne. 0) goto 1237
         read(ifile,*,iostat=ios) dtemp
         if (ios .ne. 0) goto 1237

c     first is number of modes

         read(ifile,*,iostat=ios) nmodes
         if (ios .ne. 0) goto 1237

         if (nmodes .lt. 0) then
            ier = 2
            write(*,*) 'ERROR READING FILE: '
            write(*,*) 'NMODES negative for component', i
            return
         endif

         iparstemp = iparstemp+1
         iparstemp = iparstemp+1
         iparstemp = iparstemp+1

c     then the modes (nmodes should be >= 0
c     and there should be a zero-th mode)
c     coefficients are complex, in columns of
c     two real numbers. First the x coeffs
c     are listed, then the y coeffs

c     x coeffs

         do j = 0,nmodes
            
c     real, imaginary parts of coefficients

            read(ifile,*,iostat=ios) dtemp, dtemp2
            if (ios .ne. 0) goto 1237
            iparstemp = iparstemp+1
            iparstemp = iparstemp+1         
         enddo

c     y coeffs

         do j = 0,nmodes
            
c     real, imaginary parts of coefficients

            read(ifile,*,iostat=ios) dtemp, dtemp2
            if (ios .ne. 0) goto 1237
            iparstemp = iparstemp+1
            iparstemp = iparstemp+1         
         enddo

      enddo
         
      lpars = iparstemp-1

      return

 1237 continue

      ier = 3
      write(*,*) 'ERROR READING FILE:'
      write(*,*) 'REACHED END OR BAD READ'
      write(*,*) 'CHECK FILE FORMAT'
      return

      end
