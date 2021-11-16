  
c 
c 
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c        this is the end of the debugging code and the beginning of the
c        hankel function code proper.
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
c 
        subroutine hanks104(z,hanks,n,ifexpon)
        implicit real *8 (a-h,o-z)
        save
        complex *16 z,hanks(1)
c 
c       This subroutine evaluates the first n+1 Hankel functions of the
c       argument z. The user also has the option of evaluating the
c       functions H_m(z) scaled by the (complex) coefficient e^{-i \cdot z}.
c       This option is provided via the parameter ifexpon (see below)
c 
c 
c                      input parameters:
c 
c  z - the complex number for which the hankel functions
c        H_0, H_1 are to be evaluated
c  n - the highest order of any Hankel function to be evaluated
c  ifexpon - the integer parameter telling the subroutine whether
c        to calculate the actual values of the hankel functions,
c        or the values of Hankel functions scaled by e^{-i \cdot z}.
c        Permitted values: 0 and 1.
c    ifexpon = 1 will cause the subroutine to evaluate the Hankel functions
c        honestly
c    ifexpon = 0 will cause the subroutine to scale the Hankel functions
c        by e^{-i \cdot z}.
c 
c                      output parameters:
c 
c  hanks - the first n+1 Hankel functions of the (complex) argument z.
c        Please note that hanks(1) is the Hankel function of order 0,
c        hanks(2) is the Hankel function of order 1, ..., hanks(n+1)
c        is the Hankel function of order n
c 
c       . . . evaluate the functions h0,h1
c 
        call hank103(z,hanks(1),hanks(2),ifexpon)
c 
c 
c       conduct recursion
c 
        do 1200 i1=2,n
c 
        i=i1-1
c 
        hanks(i1+1)=(2*i)/z*hanks(i1)-hanks(i1-1)
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine hanks103(z,hanks,n,ifexpon)
        implicit real *8 (a-h,o-z)
        save
        complex *16 z,hanks(1),cd,cdd
c 
c       This subroutine evaluates the first n+1 Hankel functions of the
c       argument z. The user also has the option of evaluating the
c       functions H_m(z) scaled by the (complex) coefficient e^{-i \cdot z}.
c       This option is provided via the parameter ifexpon (see below)
c 
c 
c                      input parameters:
c 
c  z - the complex number for which the hankel functions
c        H_0, H_1 are to be evaluated
c  n - the highest order of any Hankel function to be evaluated
c  ifexpon - the integer parameter telling the subroutine whether
c        to calculate the actual values of the hankel functions,
c        or the values of Hankel functions scaled by e^{-i \cdot z}.
c        Permitted values: 0 and 1.
c    ifexpon = 1 will cause the subroutine to evaluate the Hankel functions
c        honestly
c    ifexpon = 0 will cause the subroutine to scale the Hankel functions
c        by e^{-i \cdot z}.
c 
c                      output parameters:
c 
c  hanks - the first n+1 Hankel functions of the (complex) argument z.
c        Please note that hanks(1) is the Hankel function of order 0,
c        hanks(2) is the Hankel function of order 1, ..., hanks(n+1)
c        is the Hankel function of order n
c 
c       . . . evaluate the functions h0,h1
c 
        call hank103(z,hanks(1),hanks(2),ifexpon)
c 
c 
c       conduct recursion
c 
        cd=2/z
        cdd=cd
        do 1200 i1=2,n
c 
        i=i1-1
c 
cccc        hanks(i1+1)=(2*i)/z*hanks(i1)-hanks(i1-1)
        hanks(i1+1)=cdd*hanks(i1)-hanks(i1-1)
c 
        cdd=cdd+cd
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine hank103(z,h0,h1,ifexpon)
        implicit real *8 (a-h,o-z)
        save
        complex *16 z,h0,h1,h0u,h0r,h1u,h1r,
     1      fj0,fj1,y0,y1,com,zu,zr,ima,ser2,ser3,z2,
     2      cclog,cd
        real *8 rea(2)
        equivalence (rea(1),com)
        data ima/(0.0d0,1.0d0)/,pi/0.31415926535897932D+01/
c 
c        this subroutine evaluates the hankel functions H_0^1, H_1^1
c        for an arbitrary user-specified complex number z. The user
c        also has the option of evaluating the functions h0, h1
c        scaled by the (complex) coefficient e^{-i \cdot z}. This
c        subroutine is a modification of the subroutine hank102
c        (see), different from the latter by having the parameter
c        ifexpon. Please note that the subroutine hank102 is in
c        turn a slightly accelerated version of the old hank101
c        (see). The principal claim to fame of all three is that
c        they are valid on the whole  complex plane, and are
c        reasonably accurate (14-digit relative accuracy) and
c        reasonably fast. Also, please note that all three have not
c        been carefully tested in the third quadrant (both x and y
c        negative); some sort of numerical trouble is possible
c        (though has not been observed) for LARGE z in the third
c        quadrant.
c 
c                      input parameters:
c 
c  z - the complex number for which the hankel functions
c        H_0, H_1 are to be evaluated
c  ifexpon - the integer parameter telling the subroutine whether
c        to calculate the actual values of the hankel functions,
c        or the values of Hankel functions scaled by e^{-i \cdot z}.
c        Permitted values: 0 and 1.
c    ifexpon = 1 will cause the subroutine to evaluate the Hankel functions
c        honestly
c    ifexpon = 0 will cause the subroutine to scale the Hankel functions
c        by e^{-i \cdot z}.
c 
c                      output parameters:
c 
c  h0, h1 - the said Hankel functions
c 
c 
c        . . . if z in the upper half-plane - act accordingly
c 
        com=z
        if(rea(2) .lt. 0) goto 1400
        call hank103u(z,ier,h0,h1,ifexpon)
        return
 1400 continue
c 
c       if z is in the right lower quadrant - act accordingly
c 
        if(rea(1) .lt. 0) goto 2000
        call hank103r(z,ier,h0,h1,ifexpon)
        return
 2000 continue
c 
c       z is in the left lower quadrant. compute
c       h0, h1 at the points zu, zr obtained from z by reflection
c       in the x and y axis, respectively
c 
        zu=dconjg(z)
        zr=-zu
c 
        call hank103u(zu,ier,h0u,h1u,ifexpon)
        call hank103r(zr,ier,h0r,h1r,ifexpon)
  
        if(ifexpon .eq. 1) goto 3000
  
        com=zu
        subt=abs(rea(2))
  
        cd=exp(ima*zu-subt)
        h0u=h0u*cd
        h1u=h1u*cd
  
        cd=exp(ima*zr-subt)
        h0r=h0r*cd
        h1r=h1r*cd
 3000 continue
c 
c       compute the functions j0, j1, y0, y1
c       at the point zr
c 
        half=1
        half=half/2
        y0=(h0u+h0r)*half/ima
        fj0=-(h0u-h0r)*half
c 
        y1=-(h1u-h1r)*half/ima
        fj1=(h1u+h1r)*half
c 
c        finally, compute h0, h1
c 
c       . . . calculate ser2, ser3
c 
         z2=-dconjg(z)
         cclog=cdlog(z2)
         ser2=y0-fj0*2/pi*cclog
         ser3=y1-fj1*2/pi*cclog
c 
c       reflect all of these in the imaginary axis
c 
        fj0=dconjg(fj0)
        fj1=-dconjg(fj1)
c 
        ser2=dconjg(ser2)
        ser3=-dconjg(ser3)
c 
c       reconstitute y0, y1
c 
        cclog=cdlog(z)
        y0=ser2+fj0*2/pi*cclog
        y1=ser3+fj1*2/pi*cclog
c 
        h0=fj0+ima*y0
        h1=fj1+ima*y1
  
        if(ifexpon .eq. 1) return
  
        cd=exp(-ima*z+subt)
        h0=h0*cd
        h1=h1*cd
  
        return
        end
c 
c 
c 
c 
c 
        subroutine hank103u(z,ier,h0,h1,ifexpon)
        implicit real *8 (a-h,o-z)
        save
        complex *16 z,com,ima,cd,h0,h1,ccex,zzz9
        dimension rea(2)
        real *8 c0p1(34),c0p1b(36),buf01(2)
        equivalence (c0p1(34),buf01(1)),
     1      (c0p1b(1),buf01(2)),(rea(1),com)
        real *8 c1p1(34),c1p1b(36),buf11(2)
        equivalence (c1p1(34),buf11(1)),
     1      (c1p1b(1),buf11(2))
        real *8 c0p2(34),c0p2b(28),buf02(2)
        equivalence (c0p2(34),buf02(1)),
     1      (c0p2b(1),buf02(2))
        real *8 c1p2(34),c1p2b(28),buf12(2)
        equivalence (c1p2(34),buf12(1)),
     1      (c1p2b(1),buf12(2))
        data ima/(0.0d0,1.0d0)/
c 
c        this subroutine evaluates the hankel functions H_0^1, H_1^1
c        for a user-specified complex number z in the upper half-plane.
c        it is reasonably accurate (14-digit relative accuracy)
c        and reasonably fast.
c 
c 
c                      input parameters:
c 
c  z - the complex number for which the hankel functions
c        H_0, H_1 are to be evaluated
c 
c                      output parameters:
c 
c  ier - error return code.
c         ier=0 means successful conclusion
c         ier=4 means that z is not in the upper half-plane
c  h0, h1 - the said Hankel functions
c 
        data c0p1/
     1     -.6619836118357782D-12,  -.6619836118612709D-12,
     2     -.7307514264754200D-21,  0.3928160926261892D-10,
     3     0.5712712520172854D-09,  -.5712712519967086D-09,
     4     -.1083820384008718D-07,  -.1894529309455499D-18,
     5     0.7528123700585197D-07,  0.7528123700841491D-07,
     6     0.1356544045548053D-16,  -.8147940452202855D-06,
     7     -.3568198575016769D-05,  0.3568198574899888D-05,
     8     0.2592083111345422D-04,  0.4209074870019400D-15,
     9     -.7935843289157352D-04,  -.7935843289415642D-04,
     a     -.6848330800445365D-14,  0.4136028298630129D-03,
     1     0.9210433149997867D-03,  -.9210433149680665D-03,
     2     -.3495306809056563D-02,  -.6469844672213905D-13,
     3     0.5573890502766937D-02,  0.5573890503000873D-02,
     4     0.3767341857978150D-12,  -.1439178509436339D-01,
     5     -.1342403524448708D-01,  0.1342403524340215D-01,
     6     0.8733016209933828D-02,  0.1400653553627576D-11,
     7     0.2987361261932706D-01,  0.2987361261607835D-01/
        data c0p1b/
     8     -.3388096836339433D-11,  -.1690673895793793D+00,
     9     0.2838366762606121D+00,  -.2838366762542546D+00,
     a     0.7045107746587499D+00,  -.5363893133864181D-11,
     1     -.7788044738211666D+00,  -.7788044738130360D+00,
     2     0.5524779104964783D-11,  0.1146003459721775D+01,
     3     0.6930697486173089D+00,  -.6930697486240221D+00,
     4     -.7218270272305891D+00,  0.3633022466839301D-11,
     5     0.3280924142354455D+00,  0.3280924142319602D+00,
     6     -.1472323059106612D-11,  -.2608421334424268D+00,
     7     -.9031397649230536D-01,  0.9031397649339185D-01,
     8     0.5401342784296321D-01,  -.3464095071668884D-12,
     9     -.1377057052946721D-01,  -.1377057052927901D-01,
     a     0.4273263742980154D-13,  0.5877224130705015D-02,
     1     0.1022508471962664D-02,  -.1022508471978459D-02,
     2     -.2789107903871137D-03,  0.2283984571396129D-14,
     3     0.2799719727019427D-04,  0.2799719726970900D-04,
     4     -.3371218242141487D-16,  -.3682310515545645D-05,
     5     -.1191412910090512D-06,  0.1191412910113518D-06/
c 
        data c1p1/
  
     1     0.4428361927253983D-12,  -.4428361927153559D-12,
     2     -.2575693161635231D-10,  -.2878656317479645D-21,
     3     0.3658696304107867D-09,  0.3658696304188925D-09,
     4     0.7463138750413651D-19,  -.6748894854135266D-08,
     5     -.4530098210372099D-07,  0.4530098210271137D-07,
     6     0.4698787882823243D-06,  0.5343848349451927D-17,
     7     -.1948662942158171D-05,  -.1948662942204214D-05,
     8     -.1658085463182409D-15,  0.1316906100496570D-04,
     9     0.3645368564036497D-04,  -.3645368563934748D-04,
     a     -.1633458547818390D-03,  -.2697770638600506D-14,
     1     0.2816784976551660D-03,  0.2816784976676616D-03,
     2     0.2548673351180060D-13,  -.6106478245116582D-03,
     3     0.2054057459296899D-03,  -.2054057460218446D-03,
     4     -.6254962367291260D-02,  0.1484073406594994D-12,
     5     0.1952900562500057D-01,  0.1952900562457318D-01,
     6     -.5517611343746895D-12,  -.8528074392467523D-01,
     7     -.1495138141086974D+00,  0.1495138141099772D+00/
c 
        data c1p1b/
     8     0.4394907314508377D+00,  -.1334677126491326D-11,
     9     -.1113740586940341D+01,  -.1113740586937837D+01,
     a     0.2113005088866033D-11,  0.1170212831401968D+01,
     1     0.1262152242318805D+01,  -.1262152242322008D+01,
     2     -.1557810619605511D+01,  0.2176383208521897D-11,
     3     0.8560741701626648D+00,  0.8560741701600203D+00,
     4     -.1431161194996653D-11,  -.8386735092525187D+00,
     5     -.3651819176599290D+00,  0.3651819176613019D+00,
     6     0.2811692367666517D+00,  -.5799941348040361D-12,
     7     -.9494630182937280D-01,  -.9494630182894480D-01,
     8     0.1364615527772751D-12,  0.5564896498129176D-01,
     9     0.1395239688792536D-01,  -.1395239688799950D-01,
     a     -.5871314703753967D-02,  0.1683372473682212D-13,
     1     0.1009157100083457D-02,  0.1009157100077235D-02,
     2     -.8997331160162008D-15,  -.2723724213360371D-03,
     3     -.2708696587599713D-04,  0.2708696587618830D-04,
     4     0.3533092798326666D-05,  -.1328028586935163D-16,
     5     -.1134616446885126D-06,  -.1134616446876064D-06/
c 
        data c0p2/
     1     0.5641895835516786D+00,  -.5641895835516010D+00,
     2     -.3902447089770041D-09,  -.3334441074447365D-11,
     3     -.7052368835911731D-01,  -.7052368821797083D-01,
     4     0.1957299315085370D-08,  -.3126801711815631D-06,
     5     -.3967331737107949D-01,  0.3967327747706934D-01,
     6     0.6902866639752817D-04,  0.3178420816292497D-06,
     7     0.4080457166061280D-01,  0.4080045784614144D-01,
     8     -.2218731025620065D-04,  0.6518438331871517D-02,
     9     0.9798339748600499D-01,  -.9778028374972253D-01,
     a     -.3151825524811773D+00,  -.7995603166188139D-03,
     1     0.1111323666639636D+01,  0.1116791178994330D+01,
     2     0.1635711249533488D-01,  -.8527067497983841D+01,
     3     -.2595553689471247D+02,  0.2586942834408207D+02,
     4     0.1345583522428299D+03,  0.2002017907999571D+00,
     5     -.3086364384881525D+03,  -.3094609382885628D+03,
     6     -.1505974589617013D+01,  0.1250150715797207D+04,
     7     0.2205210257679573D+04,  -.2200328091885836D+04/
        data c0p2b/
     8     -.6724941072552172D+04,  -.7018887749450317D+01,
     9     0.8873498980910335D+04,  0.8891369384353965D+04,
     a     0.2008805099643591D+02,  -.2030681426035686D+05,
     1     -.2010017782384992D+05,  0.2006046282661137D+05,
     2     0.3427941581102808D+05,  0.3432892927181724D+02,
     3     -.2511417407338804D+05,  -.2516567363193558D+05,
     4     -.3318253740485142D+02,  0.3143940826027085D+05,
     5     0.1658466564673543D+05,  -.1654843151976437D+05,
     6     -.1446345041326510D+05,  -.1645433213663233D+02,
     7     0.5094709396573681D+04,  0.5106816671258367D+04,
     8     0.3470692471612145D+01,  -.2797902324245621D+04,
     9     -.5615581955514127D+03,  0.5601021281020627D+03,
     a     0.1463856702925587D+03,  0.1990076422327786D+00,
     1     -.9334741618922085D+01,  -.9361368967669095D+01/
c 
        data c1p2/
     1     -.5641895835446003D+00,  -.5641895835437973D+00,
     2     0.3473016376419171D-10,  -.3710264617214559D-09,
     3     0.2115710836381847D+00,  -.2115710851180242D+00,
     4     0.3132928887334847D-06,  0.2064187785625558D-07,
     5     -.6611954881267806D-01,  -.6611997176900310D-01,
     6     -.3386004893181560D-05,  0.7146557892862998D-04,
     7     -.5728505088320786D-01,  0.5732906930408979D-01,
     8     -.6884187195973806D-02,  -.2383737409286457D-03,
     9     0.1170452203794729D+00,  0.1192356405185651D+00,
     a     0.8652871239920498D-02,  -.3366165876561572D+00,
     1     -.1203989383538728D+01,  0.1144625888281483D+01,
     2     0.9153684260534125D+01,  0.1781426600949249D+00,
     3     -.2740411284066946D+02,  -.2834461441294877D+02,
     4     -.2192611071606340D+01,  0.1445470231392735D+03,
     5     0.3361116314072906D+03,  -.3270584743216529D+03,
     6     -.1339254798224146D+04,  -.1657618537130453D+02,
     7     0.2327097844591252D+04,  0.2380960024514808D+04/
        data c1p2b/
     8     0.7760611776965994D+02,  -.7162513471480693D+04,
     9     -.9520608696419367D+04,  0.9322604506839242D+04,
     a     0.2144033447577134D+05,  0.2230232555182369D+03,
     1     -.2087584364240919D+05,  -.2131762020653283D+05,
     2     -.3825699231499171D+03,  0.3582976792594737D+05,
     3     0.2642632405857713D+05,  -.2585137938787267D+05,
     4     -.3251446505037506D+05,  -.3710875194432116D+03,
     5     0.1683805377643986D+05,  0.1724393921722052D+05,
     6     0.1846128226280221D+03,  -.1479735877145448D+05,
     7     -.5258288893282565D+04,  0.5122237462705988D+04,
     8     0.2831540486197358D+04,  0.3905972651440027D+02,
     9     -.5562781548969544D+03,  -.5726891190727206D+03,
     a     -.2246192560136119D+01,  0.1465347141877978D+03,
     1     0.9456733342595993D+01,  -.9155767836700837D+01/
c 
c        if the user-specified z is in the lower half-plane
c        - bomb out
c 
        ier=0
        com=z
        if(rea(2) .ge. 0) goto 1200
        ier=4
        return
 1200 continue
c 
        done=1
        thresh1=1**2
        thresh2=3.7**2
        thresh3=20**2
c 
c       check if if the user-specified z is in one of the
c       intermediate regimes
c 
        d=z*dconjg(z)
        if( (d .lt. thresh1) .or. (d .gt. thresh3) ) goto 3000
c 
c        the user-specified z is in one of the intermediate regimes.
c        act accordingly
c 
c 
        if(d .gt. thresh2) goto 2000
c 
c       z is in the first intermediate regime: its absolute value is
c       between 1 and 3.7. act accordingly
c 
c       . . . evaluate the expansion
c 
        cd=done/cdsqrt(z)
c 
        ccex=cd
        if(ifexpon .eq. 1) ccex=ccex*cdexp(ima*z)
c 
        zzz9=z**9
        m=35
        call hank103p(c0p1,m,cd,h0)
        h0=h0*ccex * zzz9
c 
        call hank103p(c1p1,m,cd,h1)
        h1=h1*ccex * zzz9
        return
 2000 continue
c 
c       z is in the second intermediate regime: its absolute value is
c       between 3.7 and 20. act accordingly.
c 
        cd=done/cdsqrt(z)
c 
        ccex=cd
        if(ifexpon .eq. 1) ccex=ccex*cdexp(ima*z)
  
        m=31
        call hank103p(c0p2,m,cd,h0)
        h0=h0*ccex
c 
        m=31
        call hank103p(c1p2,m,cd,h1)
        h1=h1*ccex
        return
 3000 continue
c 
c        z is either in the local regime or the asymptotic one.
c        if it is in the local regime - act accordingly.
c 
        if(d .gt. 50.d0) goto 4000
        call hank103l(z,h0,h1,ifexpon)
        return
c 
c        z is in the asymptotic regime. act accordingly.
c 
 4000 continue
        call hank103a(z,h0,h1,ifexpon)
        return
        end
c 
c 
c 
c 
        subroutine hank103p(p,m,z,f)
        implicit real *8 (a-h,o-z)
        save
        complex *16 p(1),z,f
c 
c       evaluate a polynomial at a point
c 
        f=p(m)
        do 1200 i=m-1,1,-1
        f=f*z+p(i)
 1200 continue
        return
        end
  
  
  
  
c 
c 
c 
c 
c 
        subroutine hank103a(z,h0,h1,ifexpon)
        implicit real *8 (a-h,o-z)
        save
        dimension p(18),q(18),p1(18),q1(18),rea(2)
        complex *16 z,zinv,pp,qq,ima,h0,h1,pp1,qq1,
     1      com,cccexp,cdd,cdumb,zinv22
        equivalence (rea(1),com)
        data ima/(0.0d0,1.0d0)/,pi/0.31415926535897932D+01/,
     1      done/1.0d0/,cdumb/
     2      (0.70710678118654757D+00,-.70710678118654746D+00)/
c 
         data p/
     1     0.1000000000000000D+01,  -.7031250000000000D-01,
     2     0.1121520996093750D+00,  -.5725014209747314D+00,
     3     0.6074042001273483D+01,  -.1100171402692467D+03,
     4     0.3038090510922384D+04,  -.1188384262567833D+06,
     5     0.6252951493434797D+07,  -.4259392165047669D+09,
     6     0.3646840080706556D+11,  -.3833534661393944D+13,
     7     0.4854014686852901D+15,  -.7286857349377657D+17,
     8     0.1279721941975975D+20,  -.2599382102726235D+22,
     9     0.6046711487532401D+24,  -.1597065525294211D+27/
c 
         data q/
     1     -.1250000000000000D+00,  0.7324218750000000D-01,
     2     -.2271080017089844D+00,  0.1727727502584457D+01,
     3     -.2438052969955606D+02,  0.5513358961220206D+03,
     4     -.1825775547429317D+05,  0.8328593040162893D+06,
     5     -.5006958953198893D+08,  0.3836255180230434D+10,
     6     -.3649010818849834D+12,  0.4218971570284096D+14,
     7     -.5827244631566907D+16,  0.9476288099260110D+18,
     8     -.1792162323051699D+21,  0.3900121292034000D+23,
     9     -.9677028801069847D+25,  0.2715581773544907D+28/
  
         data p1/
     1     0.1000000000000000D+01,  0.1171875000000000D+00,
     2     -.1441955566406250D+00,  0.6765925884246826D+00,
     3     -.6883914268109947D+01,  0.1215978918765359D+03,
     4     -.3302272294480852D+04,  0.1276412726461746D+06,
     5     -.6656367718817687D+07,  0.4502786003050393D+09,
     6     -.3833857520742789D+11,  0.4011838599133198D+13,
     7     -.5060568503314726D+15,  0.7572616461117957D+17,
     8     -.1326257285320556D+20,  0.2687496750276277D+22,
     9     -.6238670582374700D+24,  0.1644739123064188D+27/
c 
         data q1/
     1     0.3750000000000000D+00,  -.1025390625000000D+00,
     2     0.2775764465332031D+00,  -.1993531733751297D+01,
     3     0.2724882731126854D+02,  -.6038440767050702D+03,
     4     0.1971837591223663D+05,  -.8902978767070679D+06,
     5     0.5310411010968522D+08,  -.4043620325107754D+10,
     6     0.3827011346598606D+12,  -.4406481417852279D+14,
     7     0.6065091351222699D+16,  -.9833883876590680D+18,
     8     0.1855045211579829D+21,  -.4027994121281017D+23,
     9     0.9974783533410457D+25,  -.2794294288720121D+28/
c 
c        evaluate the asymptotic expansion for h0,h1 at
c        the user-supplied point z, provided it is not
c        in the fourth quadrant
c 
        m=10
        zinv=done/z
c 
        pp=p(m)
        pp1=p1(m)
        zinv22=zinv**2
c 
        qq=q(m)
        qq1=q1(m)
c 
        do 1600 i=m-1,1,-1
  
        pp=pp* zinv22+p(i)
        pp1=pp1* zinv22+p1(i)
  
        qq=qq* zinv22+q(i)
        qq1=qq1* zinv22+q1(i)
 1600 continue
c 
        qq=qq*zinv
        qq1=qq1*zinv
c 
        cccexp=1
        if(ifexpon .eq. 1) cccexp=cdexp(ima*z)
c 
        cdd=cdsqrt(2/pi*zinv)
c 
        h0=pp+ima*qq
        h0=cdd*cdumb*cccexp * h0
c 
        h1=pp1+ima*qq1
        h1=-cdd*cccexp*cdumb* h1*ima
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine hank103l(z,h0,h1,ifexpon)
        implicit real *8 (a-h,o-z)
        save
        dimension cj0(16),cj1(16),ser2(16),ser2der(16)
        complex *16 z,fj0,fj1,y0,y1,h0,h1,z2,cd,ima,cdddlog
c 
        data gamma/0.5772156649015328606d+00/
        data ima/(0.0d0,1.0d0)/,pi/0.31415926535897932D+01/,
     1      two/2.0d0/
c 
c        this subroutine evaluates the hankel functions H_0^1, H_1^1
c        for a user-specified complex number z in the local regime,
c        i. e. for cdabs(z) < 1 in the upper half-plane,
c        and for cdabs(z) < 4 in the lower half-plane,
c        it is reasonably accurate (14-digit relative accuracy) and
c        reasonably fast.
c 
c                      input parameters:
c 
c  z - the complex number for which the hankel functions
c        H_0, H_1 are to be evaluated
c 
c                      output parameters:
c 
c  h0, h1 - the said Hankel functions
c 
        data cj0/
     1     0.1000000000000000D+01,  -.2500000000000000D+00,
     2     0.1562500000000000D-01,  -.4340277777777778D-03,
     3     0.6781684027777778D-05,  -.6781684027777778D-07,
     4     0.4709502797067901D-09,  -.2402807549524439D-11,
     5     0.9385966990329841D-14,  -.2896903392077112D-16,
     6     0.7242258480192779D-19,  -.1496334396734045D-21,
     7     0.2597802772107717D-24,  -.3842903509035085D-27,
     8     0.4901662639075363D-30,  -.5446291821194848D-33/
        data cj1/
     1     -.5000000000000000D+00,  0.6250000000000000D-01,
     2     -.2604166666666667D-02,  0.5425347222222222D-04,
     3     -.6781684027777778D-06,  0.5651403356481481D-08,
     4     -.3363930569334215D-10,  0.1501754718452775D-12,
     5     -.5214426105738801D-15,  0.1448451696038556D-17,
     6     -.3291935672814899D-20,  0.6234726653058522D-23,
     7     -.9991549123491221D-26,  0.1372465538941102D-28,
     8     -.1633887546358454D-31,  0.1701966194123390D-34/
        data ser2/
     1     0.2500000000000000D+00,  -.2343750000000000D-01,
     2     0.7957175925925926D-03,  -.1412850839120370D-04,
     3     0.1548484519675926D-06,  -.1153828185281636D-08,
     4     0.6230136717695511D-11,  -.2550971742728932D-13,
     5     0.8195247730999099D-16,  -.2121234517551702D-18,
     6     0.4518746345057852D-21,  -.8061529302289970D-24,
     7     0.1222094716680443D-26,  -.1593806157473552D-29,
     8     0.1807204342667468D-32,  -.1798089518115172D-35/
        data ser2der/
     1     0.5000000000000000D+00,  -.9375000000000000D-01,
     2     0.4774305555555556D-02,  -.1130280671296296D-03,
     3     0.1548484519675926D-05,  -.1384593822337963D-07,
     4     0.8722191404773715D-10,  -.4081554788366291D-12,
     5     0.1475144591579838D-14,  -.4242469035103405D-17,
     6     0.9941241959127275D-20,  -.1934767032549593D-22,
     7     0.3177446263369152D-25,  -.4462657240925946D-28,
     8     0.5421613028002404D-31,  -.5753886457968550D-34/
c 
c        evaluate j0, j1
c 
        m=16
        fj0=0
        fj1=0
        y0=0
        y1=0
        z2=z**2
        cd=1
c 
        do 1800 i=1,m
        fj0=fj0+cj0(i)*cd
        fj1=fj1+cj1(i)*cd
        y1=y1+ser2der(i)*cd
        cd=cd*z2
        y0=y0+ser2(i)*cd
 1800 continue
        fj1=-fj1*z
c 
        cdddlog=cdlog(z/two)+gamma
        y0=cdddlog*fj0+y0
        y0=two/pi*y0
c 
        y1=y1*z
c 
        y1=-cdddlog*fj1+fj0/z+y1
        y1=-y1*two/pi
c 
        h0=fj0+ima*y0
        h1=fj1+ima*y1
c 
        if(ifexpon .eq. 1) return
c 
        cd=exp(-ima*z)
        h0=h0*cd
        h1=h1*cd
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine hank103r(z,ier,h0,h1,ifexpon)
        implicit real *8 (a-h,o-z)
        save
        complex *16 z,com,ima,cd,h0,h1,cccexp,cdd,zz18
        dimension rea(2)
        real *8 c0p1(34),c0p1b(36),buf01(2)
        equivalence (c0p1(34),buf01(1)),
     1      (c0p1b(1),buf01(2)),(rea(1),com)
        real *8 c1p1(34),c1p1b(36),buf11(2)
        equivalence (c1p1(34),buf11(1)),
     1      (c1p1b(1),buf11(2))
        real *8 c0p2(34),c0p2b(20),buf02(2)
        equivalence (c0p2(34),buf02(1)),
     1      (c0p2b(1),buf02(2))
        real *8 c1p2(34),c1p2b(28),buf12(2)
        equivalence (c1p2(34),buf12(1)),
     1      (c1p2b(1),buf12(2))
        data ima/(0.0d0,1.0d0)/
c 
c        this subroutine evaluates the hankel functions H_0^1, H_1^1
c        for a user-specified complex number z in the right lower
c        quadrant. it is reasonably accurate (14-digit relative
c        accuracy) and reasonably fast.
c 
c 
c                      input parameters:
c 
c  z - the complex number for which the hankel functions
c        H_0, H_1 are to be evaluated
c 
c                      output parameters:
c 
c  ier - error return code.
c         ier=0 means successful conclusion
c         ier=4 means that z is not in the right lower quadrant
c  h0, h1 - the said Hankel functions
c 
        data c0p1/
     1     -.4268441995428495D-23,  0.4374027848105921D-23,
     2     0.9876152216238049D-23,  -.1065264808278614D-20,
     3     0.6240598085551175D-19,  0.6658529985490110D-19,
     4     -.5107210870050163D-17,  -.2931746613593983D-18,
     5     0.1611018217758854D-15,  -.1359809022054077D-15,
     6     -.7718746693707326D-15,  0.6759496139812828D-14,
     7     -.1067620915195442D-12,  -.1434699000145826D-12,
     8     0.3868453040754264D-11,  0.7061853392585180D-12,
     9     -.6220133527871203D-10,  0.3957226744337817D-10,
     a     0.3080863675628417D-09,  -.1154618431281900D-08,
     1     0.7793319486868695D-08,  0.1502570745460228D-07,
     2     -.1978090852638430D-06,  -.7396691873499030D-07,
     3     0.2175857247417038D-05,  -.8473534855334919D-06,
     4     -.1053381327609720D-04,  0.2042555121261223D-04,
     5     -.4812568848956982D-04,  -.1961519090873697D-03,
     6     0.1291714391689374D-02,  0.9234422384950050D-03,
     7     -.1113890671502769D-01,  0.9053687375483149D-03/
        data c0p1b/
     8     0.5030666896877862D-01,  -.4923119348218356D-01,
     9     0.5202355973926321D+00,  -.1705244841954454D+00,
     a     -.1134990486611273D+01,  -.1747542851820576D+01,
     1     0.8308174484970718D+01,  0.2952358687641577D+01,
     2     -.3286074510100263D+02,  0.1126542966971545D+02,
     3     0.6576015458463394D+02,  -.1006116996293757D+03,
     4     0.3216834899377392D+02,  0.3614005342307463D+03,
     5     -.6653878500833375D+03,  -.6883582242804924D+03,
     6     0.2193362007156572D+04,  0.2423724600546293D+03,
     7     -.3665925878308203D+04,  0.2474933189642588D+04,
     8     0.1987663383445796D+04,  -.7382586600895061D+04,
     9     0.4991253411017503D+04,  0.1008505017740918D+05,
     a     -.1285284928905621D+05,  -.5153674821668470D+04,
     1     0.1301656757246985D+05,  -.4821250366504323D+04,
     2     -.4982112643422311D+04,  0.9694070195648748D+04,
     3     -.1685723189234701D+04,  -.6065143678129265D+04,
     4     0.2029510635584355D+04,  0.1244402339119502D+04,
     5     -.4336682903961364D+03,  0.8923209875101459D+02/
c 
        data c1p1/
     1     -.4019450270734195D-23,  -.4819240943285824D-23,
     2     0.1087220822839791D-20,  0.1219058342725899D-21,
     3     -.7458149572694168D-19,  0.5677825613414602D-19,
     4     0.8351815799518541D-18,  -.5188585543982425D-17,
     5     0.1221075065755962D-15,  0.1789261470637227D-15,
     6     -.6829972121890858D-14,  -.1497462301804588D-14,
     7     0.1579028042950957D-12,  -.9414960303758800D-13,
     8     -.1127570848999746D-11,  0.3883137940932639D-11,
     9     -.3397569083776586D-10,  -.6779059427459179D-10,
     a     0.1149529442506273D-08,  0.4363087909873751D-09,
     1     -.1620182360840298D-07,  0.6404695607668289D-08,
     2     0.9651461037419628D-07,  -.1948572160668177D-06,
     3     0.6397881896749446D-06,  0.2318661930507743D-05,
     4     -.1983192412396578D-04,  -.1294811208715315D-04,
     5     0.2062663873080766D-03,  -.2867633324735777D-04,
     6     -.1084309075952914D-02,  0.1227880935969686D-02,
     7     0.2538406015667726D-03,  -.1153316815955356D-01/
c 
        data c1p1b/
     8     0.4520140008266983D-01,  0.5693944718258218D-01,
     9     -.9640790976658534D+00,  -.6517135574036008D+00,
     a     0.2051491829570049D+01,  -.1124151010077572D+01,
     1     -.3977380460328048D+01,  0.8200665483661009D+01,
     2     -.7950131652215817D+01,  -.3503037697046647D+02,
     3     0.9607320812492044D+02,  0.7894079689858070D+02,
     4     -.3749002890488298D+03,  -.8153831134140778D+01,
     5     0.7824282518763973D+03,  -.6035276543352174D+03,
     6     -.5004685759675768D+03,  0.2219009060854551D+04,
     7     -.2111301101664672D+04,  -.4035632271617418D+04,
     8     0.7319737262526823D+04,  0.2878734389521922D+04,
     9     -.1087404934318719D+05,  0.3945740567322783D+04,
     a     0.6727823761148537D+04,  -.1253555346597302D+05,
     1     0.3440468371829973D+04,  0.1383240926370073D+05,
     2     -.9324927373036743D+04,  -.6181580304530313D+04,
     3     0.6376198146666679D+04,  -.1033615527971958D+04,
     4     -.1497604891055181D+04,  0.1929025541588262D+04,
     5     -.4219760183545219D+02,  -.4521162915353207D+03/
c 
        data c0p2/
     1     0.5641895835569398D+00,  -.5641895835321127D+00,
     2     -.7052370223565544D-01,  -.7052369923405479D-01,
     3     -.3966909368581382D-01,  0.3966934297088857D-01,
     4     0.4130698137268744D-01,  0.4136196771522681D-01,
     5     0.6240742346896508D-01,  -.6553556513852438D-01,
     6     -.3258849904760676D-01,  -.7998036854222177D-01,
     7     -.3988006311955270D+01,  0.1327373751674479D+01,
     8     0.6121789346915312D+02,  -.9251865216627577D+02,
     9     0.4247064992018806D+03,  0.2692553333489150D+04,
     a     -.4374691601489926D+05,  -.3625248208112831D+05,
     1     0.1010975818048476D+07,  -.2859360062580096D+05,
     2     -.1138970241206912D+08,  0.1051097979526042D+08,
     3     0.2284038899211195D+08,  -.2038012515235694D+09,
     4     0.1325194353842857D+10,  0.1937443530361381D+10,
     5     -.2245999018652171D+11,  -.5998903865344352D+10,
     6     0.1793237054876609D+12,  -.8625159882306147D+11,
     7     -.5887763042735203D+12,  0.1345331284205280D+13/
c 
        data c0p2b/
     8     -.2743432269370813D+13,  -.8894942160272255D+13,
     9     0.4276463113794564D+14,  0.2665019886647781D+14,
     a     -.2280727423955498D+15,  0.3686908790553973D+14,
     1     0.5639846318168615D+15,  -.6841529051615703D+15,
     2     0.9901426799966038D+14,  0.2798406605978152D+16,
     3     -.4910062244008171D+16,  -.5126937967581805D+16,
     4     0.1387292951936756D+17,  0.1043295727224325D+16,
     5     -.1565204120687265D+17,  0.1215262806973577D+17,
     6     0.3133802397107054D+16,  -.1801394550807078D+17,
     7     0.4427598668012807D+16,  0.6923499968336864D+16/
c 
c 
        data c1p2/
     1     -.5641895835431980D+00,  -.5641895835508094D+00,
     2     0.2115710934750869D+00,  -.2115710923186134D+00,
     3     -.6611607335011594D-01,  -.6611615414079688D-01,
     4     -.5783289433408652D-01,  0.5785737744023628D-01,
     5     0.8018419623822896D-01,  0.8189816020440689D-01,
     6     0.1821045296781145D+00,  -.2179738973008740D+00,
     7     0.5544705668143094D+00,  0.2224466316444440D+01,
     8     -.8563271248520645D+02,  -.4394325758429441D+02,
     9     0.2720627547071340D+04,  -.6705390850875292D+03,
     a     -.3936221960600770D+05,  0.5791730432605451D+05,
     1     -.1976787738827811D+06,  -.1502498631245144D+07,
     2     0.2155317823990686D+08,  0.1870953796705298D+08,
     3     -.4703995711098311D+09,  0.3716595906453190D+07,
     4     0.5080557859012385D+10,  -.4534199223888966D+10,
     5     -.1064438211647413D+11,  0.8612243893745942D+11,
     6     -.5466017687785078D+12,  -.8070950386640701D+12,
     7     0.9337074941225827D+13,  0.2458379240643264D+13/
c 
        data c1p2b/
     8     -.7548692171244579D+14,  0.3751093169954336D+14,
     9     0.2460677431350039D+15,  -.5991919372881911D+15,
     a     0.1425679408434606D+16,  0.4132221939781502D+16,
     1     -.2247506469468969D+17,  -.1269771078165026D+17,
     2     0.1297336292749026D+18,  -.2802626909791308D+17,
     3     -.3467137222813017D+18,  0.4773955215582192D+18,
     4     -.2347165776580206D+18,  -.2233638097535785D+19,
     5     0.5382350866778548D+19,  0.4820328886922998D+19,
     6     -.1928978948099345D+20,  0.1575498747750907D+18,
     7     0.3049162180215152D+20,  -.2837046201123502D+20,
     8     -.5429391644354291D+19,  0.6974653380104308D+20,
     9     -.5322120857794536D+20,  -.6739879079691706D+20,
     a     0.6780343087166473D+20,  0.1053455984204666D+20,
     1     -.2218784058435737D+20,  0.1505391868530062D+20/
c 
c        if z is not in the right lower quadrant - bomb out
c 
        ier=0
        com=z
        if( (rea(1) .ge. 0) .and. (rea(2) .le. 0) ) goto 1400
        ier=4
        return
 1400 continue
c 
        done=1
        thresh1=4**2
        thresh2=8**2
        thresh3=20**2
c 
c       check if if the user-specified z is in one of the
c       intermediate regimes
c 
        d=z*dconjg(z)
        if( (d .lt. thresh1) .or. (d .gt. thresh3) ) goto 3000
c 
c        if the user-specified z is in the first intermediate regime
c        (i.e. if its absolute value is between 4 and 8), act accordingly
c 
        if(d .gt. thresh2) goto 2000
c 
        cccexp=1
        if(ifexpon .eq. 1) cccexp=cdexp(ima*z)
        cdd=done/cdsqrt(z)
        cd=done/z
        zz18=z**18
        m=35
        call hank103p(c0p1,m,cd,h0)
        h0=h0*cdd*cccexp*zz18
c 
        call hank103p(c1p1,m,cd,h1)
        h1=h1*cdd*cccexp*zz18
        return
 2000 continue
c 
c       z is in the second intermediate regime (i.e. its
c       absolute value is between 8 and 20). act accordingly.
c 
        cd=done/z
        cdd=sqrt(cd)
  
        cccexp=1
        if(ifexpon .eq. 1) cccexp=cdexp(ima*z)
  
        m=27
c 
        call hank103p(c0p2,m,cd,h0)
        h0=h0*cccexp*cdd
c 
        m=31
        call hank103p(c1p2,m,cd,h1)
        h1=h1*cccexp*cdd
        return
 3000 continue
c 
c 
c        z is either in the local regime or the asymptotic one.
c        if it is in the local regime - act accordingly.
c 
        if(d .gt. 50.d0) goto 4000
        call hank103l(z,h0,h1,ifexpon)
        return
c 
c        z is in the asymptotic regime. act accordingly.
c 
 4000 continue
        call hank103a(z,h0,h1,ifexpon)
        return
        end
  
  
  
  
