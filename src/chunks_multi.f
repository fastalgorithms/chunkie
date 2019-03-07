
      subroutine multichunk_init(wgeos,lwgeos,ncompmax,ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This subroutine initializes a multichunk array.
c
c     INPUT
c
c     lwgeos - the length of the wgeos array
c     ncompmax - the maximum number of components to be stored in the
c                array. This affects how items are stored in the array
c                so make this number as small as possible for your
c                application.
c     
c     OUTPUT
c
c     wgeos - the multichunk array, initialized with zero components
c     ier - error flag. if ier = 0, normal execution. if ier = 1,
c           lwgeos is too short (far too short, really)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real *8 wgeos(*)
      integer lwgeos, ier, ncompmax
c     local variables
      integer ncomp, indeces, lused
c     
c     store multiple component geometries
c     

      ier = 0

      ncomp = 0
      indeces = 6
      lused = indeces+ncompmax-1

      if (lwgeos .lt. lused) then
         ier = 1
         return
      endif

      wgeos(1) = ncomp+.1d0
      wgeos(2) = indeces+.1d0
      wgeos(3) = lwgeos+.1d0
      wgeos(4) = ncompmax+.1d0
      wgeos(5) = lused+.1d0
      
      return
      end

      subroutine multichunk_append(wgeos,wgeo,ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This subroutine adds a component to the multichunk array.
c
c     INPUT
c
c     wgeos - real *8 array. multichunk that's at least been initialied.
c     wgeo - real *8 array. the packed version of the component to 
c            add to the multichunk.
c
c     OUTPUT
c
c     wgeos - real *8 array. multichunk with wgeo added to it.
c     ier - integer, flag. error flag. if ier = 0, normal execution.
c           if ier = 1, wgeos array is too short.
c           if ier = 2, wgeos array was only created to store a 
c                       certain number of components, which would
c                       be exceeded by appending this component.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real *8 wgeos(*), wgeo(*)
      integer ier
c     local variables
      integer ncomp, indeces, lwgeos, ncompmax
      integer ncomp1, lused1, lused
      integer lusednew, i
c     
c     store multiple component geometries
c     
      ier = 0

      ncomp = wgeos(1)
      indeces = wgeos(2)
      lwgeos = wgeos(3)
      ncompmax = wgeos(4)
      lused = wgeos(5)

      call chunkused(wgeo,lusednew)

      lused1 = lused+lusednew

      if(lwgeos .lt. lused1) then
         ier = 1
         return
      endif

      ncomp1 = ncomp+1

      if(ncompmax .lt. ncomp1) then
         ier = 2
         return
      endif

      wgeos(1) = ncomp1+.1d0
      wgeos(5) = lused1+.1d0

c     starting index of new component

      wgeos(indeces+ncomp1-1) = lused+1+0.1d0

c     copy over new component

      do i = 1,lusednew
         wgeos(i+lused) = wgeo(i)
      enddo
      
      return
      end

      subroutine multichunkgetcomp(wgeos,wgeo,lwgeo,icomp,ier)
      implicit none
      real *8 wgeos(*), wgeo(*)
      integer lwgeo, icomp, ier
c     local variables
      integer ncomp, indeces, index, lusedi, i

      ier = 0
      
      ncomp = wgeos(1)

      if (ncomp .lt. icomp) then
         ier = 2
         return
      endif

      indeces = wgeos(2)
      index = wgeos(indeces+icomp-1)

      call chunkused(wgeos(index),lusedi)

      if (lwgeo .lt. lusedi) then
         ier = 1
         return
      endif

      do i = 1,lusedi
         wgeo(i) = wgeos(index+i-1)
      enddo

      return
      end

      subroutine multichunkmergepack(wgeos,wgeo,lwgeo,lused,ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This routine takes the chunks from each of the components
c     stored in wgeos and forms chunks, ders, ders2, adjs, hs arrays
c     which contain all of the chunks and packs them into wgeo.
c
c     The chunks are stored so that the chunks of the first component
c     are first, the second component are second, etc.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real *8 wgeos(*), wgeo(*)
      integer lwgeo, ier
c     local variables
      real *8, allocatable :: chunks(:,:,:), ders(:,:,:), ders2(:,:,:)
      real *8, allocatable :: hs(:)
      integer, allocatable :: adjs(:,:)
      integer imode, k, nch, lused

      ier = 0

c     query number of chunks and order of chunks

      imode = 0
      call multichunk_merge(imode,wgeos,k,nch,chunks,adjs,ders,
     1     ders2,hs,ier)
      
      if (ier .ne. 0) return

c     allocate chunks, etc

      allocate(chunks(2,k,nch))
      allocate(ders(2,k,nch))
      allocate(ders2(2,k,nch))
      allocate(hs(nch))
      allocate(adjs(2,nch))

c     get chunks, etc

      imode = 1
      call multichunk_merge(imode,wgeos,k,nch,chunks,adjs,ders,
     1     ders2,hs,ier)

      call chunkpack(k,nch,chunks,adjs,ders,ders2,
     1      hs,wgeo,lused)

      if (lused .gt. lwgeo) then
         ier = 31
         return
      endif

      return
      end

      subroutine multichunk_merge(imode,wgeos,k,nch,chunks,adjs,ders,
     1     ders2,hs,ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This routine takes the chunks from each of the components
c     stored in wgeos and forms chunks, ders, ders2, adjs, hs arrays
c     which contain all of the chunks.
c
c     The chunks are stored so that the chunks of the first component
c     are first, the second component are second, etc.
c
c     INPUT
c
c     imode - integer, flag. if imode = 0, then the total number of
c             chunks nch and the order k are returned. if imode \neq 0
c             then the arrays chunks, adjs, etc are computed.
c
c     wgeos - real *8 array. multichunk
c
c     OUTPUT
c
c     if imode = 0
c
c     nch - integer, total number of chunks
c     k - integer, order of chunks
c     ier - if ier = 11, the chunk order is not consistent across
c           components. if ier = 0, normal execution
c
c     if imode \neq 0
c
c     nch - integer, total number of chunks
c     k - integer, order of chunks
c     ier - if ier = 11, the chunk order is not consistent across
c           components. if ier = 0, normal execution
c     chunks, adjs, ders, ders2, hs - in standard chunks format
c           but representing the chunks for all components.
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer imode
      real *8 wgeos(*)
      integer k, nch, ier, adjs(2,*)
      real *8 chunks(2,*), ders(2,*), ders2(2,*), hs(*)
c     local variables
      integer ncomp, ncompmax, indeces, lwgeos, lused
      integer ktemp, ktemp2, icomp, index, nchtemp
      integer nch1,ichunks1,iadjs1,iders1,iders12,ihs1
      integer, allocatable :: icompch(:)

      ier = 0

      call multichunkinfo(wgeos,ncomp,ncompmax,indeces,
     1     lwgeos,lused)

c     figure out order of chunks
c     if they're not all the same, abort

c     figure out total number of chunks

      nch = 0
      nchtemp = 0

      do icomp = 1,ncomp
         call multichunk_getcomp_index(wgeos,icomp,index,ier)
         if (ier .ne. 0) return
         call chunkunpack1(wgeos(index),ktemp2,nch1,ichunks1,
     1        iadjs1,iders1,iders12,ihs1)
         nchtemp = nchtemp+nch1
         if (icomp .gt. 1) then
            if (ktemp .ne. ktemp2) then
               ier = 11
               return
            endif
         else
            ktemp = ktemp2
         endif
      enddo

      k = ktemp

      if (imode .eq. 0) then
         nch = nchtemp
         return
      endif

      allocate(icompch(nchtemp))

      call multichunk_merge1(wgeos,k,nch,chunks,adjs,ders,
     1     ders2,hs,icompch,ier)
      
      return 
      end

      subroutine multichunk_merge_wicompch(imode,wgeos,k,nch,chunks,
     1     adjs,ders,ders2,hs,icompch,ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This routine takes the chunks from each of the components
c     stored in wgeos and forms chunks, ders, ders2, adjs, hs arrays
c     which contain all of the chunks.
c
c     The chunks are stored so that the chunks of the first component
c     are first, the second component are second, etc.
c
c     INPUT
c
c     imode - integer, flag. if imode = 0, then the total number of
c             chunks nch and the order k are returned. if imode \neq 0
c             then the arrays chunks, adjs, etc are computed.
c
c     wgeos - real *8 array. multichunk
c
c     OUTPUT
c
c     if imode = 0
c
c     nch - integer, total number of chunks
c     k - integer, order of chunks
c     ier - if ier = 11, the chunk order is not consistent across
c           components. if ier = 0, normal execution
c
c     if imode \neq 0
c
c     nch - integer, total number of chunks
c     k - integer, order of chunks
c     ier - if ier = 11, the chunk order is not consistent across
c           components. if ier = 0, normal execution
c     chunks, adjs, ders, ders2, hs - in standard chunks format
c            but representing the chunks for all components.
c     icompch - which component a given chunk belongs to
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer imode
      real *8 wgeos(*)
      integer k, nch, ier, adjs(2,*), icompch(*)
      real *8 chunks(2,*), ders(2,*), ders2(2,*), hs(*)
c     local variables
      integer ncomp, ncompmax, indeces, lwgeos, lused
      integer ktemp, ktemp2, icomp, index, nchtemp
      integer nch1,ichunks1,iadjs1,iders1,iders12,ihs1

      ier = 0

      call multichunkinfo(wgeos,ncomp,ncompmax,indeces,
     1     lwgeos,lused)

c     figure out order of chunks
c     if they're not all the same, abort

c     figure out total number of chunks

      nch = 0
      nchtemp = 0

      do icomp = 1,ncomp
         call multichunk_getcomp_index(wgeos,icomp,index,ier)
         if (ier .ne. 0) return
         call chunkunpack1(wgeos(index),ktemp2,nch1,ichunks1,
     1        iadjs1,iders1,iders12,ihs1)
         nchtemp = nchtemp+nch1
         if (icomp .gt. 1) then
            if (ktemp .ne. ktemp2) then
               ier = 11
               return
            endif
         else
            ktemp = ktemp2
         endif
      enddo

      k = ktemp

      if (imode .eq. 0) then
         nch = nchtemp
         return
      endif

      call multichunk_merge1(wgeos,k,nch,chunks,adjs,ders,
     1     ders2,hs,icompch,ier)
      
      return 
      end

      subroutine multichunk_merge1(wgeos,k,nch,chunks,adjs,ders,
     1     ders2,hs,icompch,ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This subroutine is called by multichunk_merge and does most
c     of the work.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real *8 wgeos(*)
      integer k, nch, ier, adjs(2,*), icompch(*)
      real *8 chunks(2,k,*), ders(2,k,*), ders2(2,k,*), hs(*)
c     local variables
      integer ncomp, ncompmax, indeces, lwgeos, lused
      integer k1, icomp, index
      integer nch1, i

      ier = 0

      call multichunkinfo(wgeos,ncomp,ncompmax,indeces,
     1     lwgeos,lused)

c     figure out order of chunks
c     if they're not all the same, abort

c     figure out total number of chunks

      nch = 0

      do icomp = 1,ncomp
         call multichunk_getcomp_index(wgeos,icomp,index,ier)
         if (ier .ne. 0) return
         call chunkunpack(wgeos(index),k1,nch1,chunks(1,1,nch+1),
     1        adjs(1,nch+1),ders(1,1,nch+1),ders2(1,1,nch+1),
     1        hs(nch+1))
         do i = 1,nch1
            adjs(1,nch+i) = adjs(1,nch+i)+nch
            adjs(2,nch+i) = adjs(2,nch+i)+nch
            icompch(nch+i) = icomp
         enddo            
         nch = nch+nch1
      enddo

      return
      end

      subroutine multichunkinfo(wgeos,ncomp,ncompmax,indeces,
     1     lwgeos,lused)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     this subroutine returns some basic info about the multichunk
c     wgeos.
c
c     INPUT
c
c     wgeos - real *8 array. multichunk
c
c     OUTPUT
c
c     ncomp - integer. number of components in multichunk
c     ncompmax - integer. the maximum number of components the 
c                multichunk was initialized with.
c     indeces - integer. the starting point for storing the 
c               indeces of the ncomp wgeo arrays within wgeos.
c     lwgeos - integer. the length of the wgeos array (provided
c              by the user when wgeos was initialized)
c     lused - integer. the amount of the wgeos array that's been used.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer ncomp, indeces, lwgeos, ncompmax, lused
      real *8 wgeos(*)

      ncomp = wgeos(1)
      indeces = wgeos(2)
      lwgeos = wgeos(3)
      ncompmax = wgeos(4)
      lused = wgeos(5)

      return
      end

      subroutine multichunknchs(wgeos,nchs,ncomp)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This subroutine returns the number of chunks for each component
c
c     INPUT 
c     
c     wgeos - real *8 array. multichunk
c     
c     OUTPUT
c
c     nchs - integer array. nch(i) is the number of chunks on the ith
c            component
c     ncomp - integer. the number of components
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      integer nchs(*), ncomp
      real *8 wgeos(*)
c     local variables
      integer ncompmax, indeces, lwgeos, lused
      integer icomp, ier, index
      integer ktemp1,nch1,ichunks1,iadjs1,iders1,iders12,ihs1

      call multichunkinfo(wgeos,ncomp,ncompmax,indeces,
     1     lwgeos,lused)

      do icomp = 1,ncomp
         call multichunk_getcomp_index(wgeos,icomp,index,ier)
         if (ier .ne. 0) return
         call chunkunpack1(wgeos(index),ktemp1,nch1,ichunks1,
     1        iadjs1,iders1,iders12,ihs1)
         nchs(icomp) = nch1
      enddo

      return
      end

      subroutine multichunk_getcomp_index(wgeos,icomp,index,ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This subroutine returns the index within wgeos such that the
c     array beginning with wgeos(index) is the packed version of
c     the icomp-th component. 
c
c     INPUT
c
c     wgeos - real*8 array. multichunk
c     icomp - integer. the component whose index is being retrieved
c     
c     OUTPUT
c
c     index - integer. the index as described above
c     ier - integer, flag. error flag. if ier = 0, normal execution.
c           if ier = 2, then there is no icomp-th component.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real *8 wgeos(*)
      integer icomp, ier, index
c     local variables
      integer ncomp, indeces
      
      ier = 0

      ncomp = wgeos(1)

      if (ncomp .lt. icomp) then
         ier = 2
         return
      endif

      indeces = wgeos(2)
      index = wgeos(indeces+icomp-1)

      return
      end
      
      subroutine multichunkfunc(eps, ifclosed, chsmall, tas, tbs,
     1     funcurve, pars, iparsindeces, k, novers, wgeos, lwgeos,
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
c     funcurve - external subroutine with the calling sequence
c     
c           funcurve(t,pars,x,y,dxdt,dydt,dxdt2,dydt2)
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
      external funcurve
c     local variables
      integer lusedi, i
      integer ncomp1, ncompmax1, lwgeos1, indeces1
      real *8, allocatable :: wgeotemp(:)

      ier = 0

      if (ncomp .gt. ncompmax) then
         ier = 1
         return
      endif

      call multichunk_init(wgeos,lwgeos,ncompmax,ier)
      if (ier .ne. 0) then
         ier = 2
         return
      endif

      lused = wgeos(5)

      allocate(wgeotemp(lwgeos))

      do i = 1,ncomp

         call chunkfunc_pack(eps(i),ifclosed(i),chsmall(i),
     1        tas(i),tbs(i),funcurve,pars(iparsindeces(i)),
     2        k,novers(i),wgeotemp,lusedi)
         if (lusedi .gt. lwgeos) then
            ier = 11
            return
         endif
         call chunksort_pack(wgeotemp)
         call multichunk_append(wgeos,wgeotemp,ier)
         if (ier .ne. 0) then
            ier = 20+i
            return
         endif
      enddo

      call multichunkinfo(wgeos,ncomp1,ncompmax1,indeces1,
     1     lwgeos1,lused)
      
      return
      end


      subroutine multichunk_reverse(wgeos,ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     this subroutine goes through and reverses the geometry for each 
c     component
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real *8 wgeos(*)
      integer ier
c     local variables
      integer ncomp, ncompmax, indeces, lwgeos, lused
      integer ktemp, ktemp2, icomp, index, nchtemp
      integer nch1,ichunks1,iadjs1,iders1,iders12,ihs1

      ier = 0

      call multichunkinfo(wgeos,ncomp,ncompmax,indeces,
     1     lwgeos,lused)

c     go through components, reversing each

      do icomp = 1,ncomp
         call multichunk_getcomp_index(wgeos,icomp,index,ier)
         if (ier .ne. 0) return
         call chunkunpack1(wgeos(index),ktemp2,nch1,ichunks1,
     1        iadjs1,iders1,iders12,ihs1)
         call chunkreverse(ktemp2,nch1,wgeos(index+ichunks1-1),
     1        wgeos(index+iadjs1-1),wgeos(index+iders1-1),
     2        wgeos(index+iders12-1),wgeos(index+ihs1-1))
      enddo

      return 
      end

      subroutine multichunk_reverse_comp(wgeos,icomp,ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     this subroutine reverses the geometry for only the icomp-th 
c     component
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      real *8 wgeos(*)
      integer ier, icomp
c     local variables
      integer ncomp, ncompmax, indeces, lwgeos, lused
      integer ktemp, ktemp2, index, nchtemp
      integer nch1,ichunks1,iadjs1,iders1,iders12,ihs1

      ier = 0

      call multichunkinfo(wgeos,ncomp,ncompmax,indeces,
     1     lwgeos,lused)

c     reversing only this component

      call multichunk_getcomp_index(wgeos,icomp,index,ier)
      if (ier .ne. 0) return
      call chunkunpack1(wgeos(index),ktemp2,nch1,ichunks1,
     1     iadjs1,iders1,iders12,ihs1)
      call chunkreverse(ktemp2,nch1,wgeos(index+ichunks1-1),
     1     wgeos(index+iadjs1-1),wgeos(index+iders1-1),
     2     wgeos(index+iders12-1),wgeos(index+ihs1-1))

      return 
      end

      
      subroutine multichunksort(wgeos,ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     this subroutine goes through and sorts the geometry for each 
c     component
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real *8 wgeos(*)
      integer ier
c     local variables
      integer ncomp, ncompmax, indeces, lwgeos, lused
      integer ktemp, ktemp2, icomp, index, nchtemp
      integer nch1,ichunks1,iadjs1,iders1,iders12,ihs1

      ier = 0

      call multichunkinfo(wgeos,ncomp,ncompmax,indeces,
     1     lwgeos,lused)

c     go through components, sorting each

      do icomp = 1,ncomp
         call multichunk_getcomp_index(wgeos,icomp,index,ier)
         if (ier .ne. 0) return
         call chunksort_pack(wgeos(index))
      enddo

      return 
      end

      subroutine multichunk_sort_comp(wgeos,icomp,ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     this subroutine sorts the geometry for only the icomp-th 
c     component
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real *8 wgeos(*)
      integer ier, icomp
c     local variables
      integer ncomp, ncompmax, indeces, lwgeos, lused
      integer ktemp, ktemp2, index, nchtemp
      integer nch1,ichunks1,iadjs1,iders1,iders12,ihs1

      ier = 0

      call multichunkinfo(wgeos,ncomp,ncompmax,indeces,
     1     lwgeos,lused)

c     sorting only this component

      call multichunk_getcomp_index(wgeos,icomp,index,ier)
      if (ier .ne. 0) return
      call chunksort_pack(wgeos(index))

      return 
      end

      
