c         This file contains four two user-callable subroutines: 
c         gammanew_eval, gammanew_eval_extend, gammanew_eval_log, 
c         gammanew_eval_extend_log. The first calculates the Gamma 
c         function of a real argument (either positive or negative) 
c         to full double precision. The second calculates the Gamma 
c         function of a real argument (either positive or negative) 
c         to (almost) full extended precision, provided the variables 
c         are re-declared to support such precision.
c
c         The subroutines gammanew_eval_log, gammanew_eval_extend_log
c         are identical to the subroutines gammanew_eval, 
c         gammanew_eval_extend, with the exception that instead of the
c         gamma function they return the said function's logarithm.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
        subroutine gammanew_eval_log_extend(x,gam)
        implicit real *8 (a-h,o-z)
        save
c
c        This subroutine evaluates the logarithm of the Gamma function 
c        of a real argument; it produces (almost) full extended 
c        precision results for both positive and negative values 
c        of the argument, provided the variables are re-declared 
c        to support such precision.
c
c        PLEASE NOTE THAT IN THE REGIMES WHERE THE GAMMA FUNCTION IS
C        NEGATIVE, THE SUBROUTINE RETURNS THE LOGARITHM OF ITS 
C        ABSOLUTE VALUE.
c
c                  Input parameters:
c
c  x - the argument of which the Gamma function is to be computed
c
c                  Output parameters:
c
c  gam - the logarithm of the Gamma function of x (how very appropriate!)
c
c
c       . . . evaluate the logarithm of the Gamma function 
c             if the argument is on the interval [1,2]
c
        if( (x .lt. 1) .or. (x.gt. 2) ) goto 2200
c
        call gammanew_eval3_log(x,gam)
        return
c
 2200 continue
c
c
c       ... if x > 2
c
        if(x .lt. 1.5) goto 3200
c
c       find x0 on the interval [1,2] with which we start,
c       and evaluate gamma at it
c
        n=x
        x0=x-n+1
        call gammanew_eval3_log(x0,gam)
c
c       march!
c
        do 2400 i=2,n
c
        gam=gam+log(i+x0-2)
 2400 continue
c
        return
c
 3200 continue
c
c       ... if  0 .leq. x .leq.  1
c
        if(x .lt. 0) goto 4200
c
        x0=x+1
        call gammanew_eval3_log(x0,gam)
        gam=gam-log(x)
        return
c
 4200 continue
c
c       ... if x < 0
c
        x1=-x
        n=x1
        x0=x1-n
        x0=1-x0+1
c
        call gammanew_eval3_log(x0,gam)
c
        xx=x0-1
        do 5200 i=1,n+2
c
        gam=gam-log(abs(xx))
        xx=xx-1
 5200 continue
c
        return
        end        
c
c
c
c
c
        subroutine gammanew_eval_log(x,gam)
        implicit real *8 (a-h,o-z)
        save
c
c        This subroutine evaluates the logarithm of the Gamma function 
c        of a real argument; it produces full double precision results 
c        for both positive and negative values of the argument.
c
C        PLEASE NOTE THAT IN THE REGIMES WHERE THE GAMMA FUNCTION IS
C        NEGATIVE, THE SUBROUTINE RETURNS THE LOGARITHM OF ITS 
C        ABSOLUTE VALUE.
C
c                  Input parameters:
c
c  x - the argument of which the Gamma function is to be computed
c
c                  Output parameters:
c
c  gam - the logarithm of the Gamma function of x (how very appropriate!)
c
c
c       . . . evaluate the logarithm of the Gamma function 
c             if the argument is on the interval [1,2]
c
        if( (x .lt. 1) .or. (x.gt. 2) ) goto 2200
c
        call gammanew_eval2_log(x,gam)
        return
c
 2200 continue
c
c
c       ... if x > 2
c
        if(x .lt. 1.5) goto 3200
c
c       find x0 on the interval [1,2] with which we start,
c       and evaluate gamma at it
c
        n=x
        x0=x-n+1
        call gammanew_eval2_log(x0,gam)
c
c       march!
c
        do 2400 i=2,n
c
        gam=gam+log(i+x0-2)
 2400 continue
c
        return
c
 3200 continue
c
c       ... if  0 .leq. x .leq.  1
c
        if(x .lt. 0) goto 4200
c
        x0=x+1
        call gammanew_eval2_log(x0,gam)
        gam=gam-log(x)
        return
c
 4200 continue
c
c       ... if x < 0
c
        x1=-x
        n=x1
        x0=x1-n
        x0=1-x0+1
c
        call prin2('x0=*',x0,1)

        call gammanew_eval2_log(x0,gam)
c
        xx=x0-1
        do 5200 i=1,n+2
c
        gam=gam-log(abs(xx))
        xx=xx-1
 5200 continue
c
        return
        end        
c
c 
c 
c 
c 
        subroutine gammanew_eval3_log(x,gam)
        implicit real *8 (a-h,o-z)
        dimension coefs(42)
        data coefs/
     1  -.120782237635245222345518445781646D+00,
     2  0.182449869892882602795118335006149D-01,
     3  0.116850275068084913677155687492027D+00,
     4  -.172665967548816665749236304400537D-01,
     5  0.366950790104801363656336640942244D-02,
     6  -.904752559027923226702063033929312D-03,
     7  0.241179440157020317746429313847828D-03,
     8  -.673640931966506793016607029366479D-04,
     9  0.193973781620149128734637856623748D-04,
     *  -.570502042708584364645567397395294D-05,
     1  0.170413630448254881730170850280867D-05,
     2  -.515095553646304672163231216170731D-06,
     3  0.157154048593298006175111322327474D-06,
     4  -.483119555232553485171697583494630D-07,
     5  0.149457513722498842477000858912645D-07,
     6  -.464831354195249220547668790930112D-08,
     7  0.145232233619810644923193813836471D-08,
     8  -.455578791522056192987778762139738D-09,
     9  0.143413197588172550414382835362441D-09,
     *  -.452865323417441015573380473066286D-10,
     1  0.143403848882077097551112487075603D-10,
     2  -.455243644934763089569345512431042D-11,
     3  0.144848973320723395282144093531243D-11,
     4  -.461834869766318701759601704380834D-12,
     5  0.147530229774684050616043905322468D-12,
     6  -.472095899923600758417333684711329D-13,
     7  0.151310521710013342466880922231033D-13,
     8  -.485686801476166232271748203132829D-14,
     9  0.156144198997796652518964422860473D-14,
     *  -.502544804552392246881299973067650D-15,
     1  0.161580172500207145883271023448381D-15,
     2  -.521114590948289341160665737783119D-16,
     3  0.171382083302691489326179629126128D-16,
     4  -.554748821197081512427887357740726D-17,
     5  0.158673415478853248945912441533215D-17,
     6  -.509788453463554343682988517338423D-18,
     7  0.267137492143822801798879951011477D-18,
     8  -.880385756441013041066119183271692D-19,
     9  -.586092705483764190954532218784502D-20,
     *  0.220231881515274593029003021427310D-20,
     1  0.645061720132042199184763314365149D-20,
     2  -.212720240682484916622372020931097D-20/
c
        save
c
        n=42
        done=1
        gam=0
c
        t=(x-1-done/2)*2 
c
        do 1200 i=1,n
c
        gam=gam+coefs(i)*t**(i-1)
cccc        gam=gam+coefs(i)*x**(i-1)
 1200 continue
c
        return
        end
c
c 
c 
c 
c 
        subroutine gammanew_eval2_log(x,gam)
        implicit real *8 (a-h,o-z)
        dimension coefs(22)
c
        data coefs/
     1  -.120782237635245221D+00,0.182449869892882600D-01,
     2  0.116850275068084696D+00,-.172665967548815972D-01,
     3  0.366950790105704649D-02,-.904752559030809163D-03,
     4  0.241179440011251918D-03,-.673640931500792339D-04,
     5  0.193973793625582700D-04,-.570502081063888534D-05,
     6  0.170413056607272862D-05,-.515093720376700298D-06,
     7  0.157171045376283400D-06,-.483173853131163447D-07,
     8  0.149136746086145614D-07,-.463806694457011235D-08,
     9  0.149074784822671166D-08,-.467852246675690637D-09,
     *  0.115448734319838611D-09,-.363558981433883658D-10,
     1  0.252470225044276880D-10,-.803435452853097061D-11/
c
        save
c
        n=22
        done=1
        gam=0
c
        t=(x-1-done/2)*2 
c
        do 1200 i=1,n
c
        gam=gam+coefs(i)*t**(i-1)
cccc        gam=gam+coefs(i)*x**(i-1)
 1200 continue
c
        return
        end
c
c
c
c
c
        subroutine gammanew_eval(x,gam)
        implicit real *8 (a-h,o-z)
        save
c
c        This subroutine evaluates the Gamma function of a real argument;
c        it produces full double precision results for both positive and 
c        negative values of the argument.
c
c                  Input parameters:
c
c  x - the argument of which the Gamma function is to be computed
c
c                  Output parameters:
c
c  gam - the Gamma function of x (how very appropriate!)
c
c
c       . . . evaluate the Gamma function if the argument is
c             on the interval [1,2]
c
        if( (x .lt. 1) .or. (x.gt. 2) ) goto 2200
c
        call gammanew_eval2(x,gam)
        return
c
 2200 continue
c
c
c       ... if the x > 2
c
        if(x .lt. 1.5) goto 3200
c
c       find x0 on the interval [1,2] with which we start,
c       and evaluate gamma at it
c
        n=x
        x0=x-n+1
        call gammanew_eval2(x0,gam)
c
c       march!
c
        do 2400 i=2,n
c
        gam=gam*(i+x0-2)
 2400 continue
c
        return
c
 3200 continue
c
c       ... if the 0 .leq. x .leq.  1
c
        if(x .lt. 0) goto 4200
c
        x0=x+1
        call gammanew_eval2(x0,gam)
        gam=gam/x
        return
c
 4200 continue
c
c       ... if x < 0
c
        x1=-x
        n=x1
        x0=x1-n
        x0=1-x0+1
c
        call gammanew_eval2(x0,gam)
c
        xx=x0-1
        do 5200 i=1,n+2
c
        gam=gam/xx
        xx=xx-1
 5200 continue
c
        return
        end        
c
c 
c 
c 
c 
        subroutine gammanew_eval_extend(x,gam)
        implicit real *8 (a-h,o-z)
        save
c
c        This subroutine evaluates the Gamma function of a real argument;
c        it produces (almost) full extended precision results for both 
c        positive and negative values of the argument, provided the 
c        variables are re-declared to support such precision.
c
c                  Input parameters:
c
c  x - the argument of which the Gamma function is to be computed
c
c                  Output parameters:
c
c  gam - the Gamma function of x (how very appropriate!)
c
c
c
c       evaluate the Gamma function if the argument is
c       on the interval [1,2]
c
        if( (x .lt. 1) .or. (x.gt. 2) ) goto 2200
c
        call gammanew_eval3(x,gam)
        return
c
 2200 continue
c
c
c       ... if the x > 2
c
        if(x .lt. 1.5) goto 3200
c
c       find x0 on the interval [1,2] with which we start,
c       and evaluate gamma at it
c
        n=x
        x0=x-n+1
        call gammanew_eval3(x0,gam)
c
c       march!
c
        do 2400 i=2,n
c
        gam=gam*(i+x0-2)
 2400 continue
c
        return
c
 3200 continue
c
c       ... if the 0 .leq. x .leq.  1
c
        if(x .lt. 0) goto 4200
c
        x0=x+1
        call gammanew_eval3(x0,gam)
        gam=gam/x
        return
c
 4200 continue
c
c       ... if x < 0
c
        x1=-x
        n=x1
        x0=x1-n
        x0=1-x0+1
c
        call gammanew_eval3(x0,gam)
c
        xx=x0-1
        do 5200 i=1,n+2
c
        gam=gam/xx
        xx=xx-1
 5200 continue
c
        return
        end        
c
c 
c 
c 
c 
        subroutine gammanew_eval3(x,gam)
        implicit real *8 (a-h,o-z)
        dimension coefs(42)
c
        data coefs/
     1  0.886226925452758013649083741670623D+00,
     2  0.161691987244425069144349421339960D-01,
     3  0.103703363422075292057509405775161D+00,
     4  -.134118505705965264609427444590074D-01,
     5  0.904033494028884643989576361780254D-02,
     6  -.242259538437044385764617259435110D-02,
     7  0.915785997143379523502915231849406D-03,
     8  -.296890121522383830093923366409561D-03,
     9  0.100928150217797671460960524521708D-03,
     *  -.336375842059855958365225717995453D-04,
     1  0.112524564378898714084403565633111D-04,
     2  -.375498601769608326351463199728310D-05,
     3  0.125284265183407932726509306293675D-05,
     4  -.417822570478480417873660426799984D-06,
     5  0.139318222887586378078540067572924D-06,
     6  -.464480612525716981485156868788579D-07,
     7  0.154844321007218425984980044056116D-07,
     8  -.516182587481302823323689251947993D-08,
     9  0.172067842623297864468229367955008D-08,
     *  -.573573438914607316487896043927398D-09,
     1  0.191193940645124209184837634463979D-09,
     2  -.637318726733661744876732845701618D-10,
     3  0.212440679157787075404568786925963D-10,
     4  -.708137779578620656802105580438353D-11,
     5  0.236046717483456631493100121338592D-11,
     6  -.786824290629355714996723594265115D-12,
     7  0.262268505662015465618725444423847D-12,
     8  -.874213989908817094617396890755866D-13,
     9  0.291499417333710813797738200069174D-13,
     *  -.971833531575606294591946449847318D-14,
     1  0.322856823583449058974923602843996D-14,
     2  -.107469659233788012193712173282348D-14,
     3  0.367876256573690494365473260918666D-15,
     4  -.123625276708537256683380049874862D-15,
     5  0.347433987117776040343629749755013D-16,
     6  -.110914923764776138469573598500169D-16,
     7  0.686389384962107642842040168449156D-17,
     8  -.245328327855757521477785550619204D-17,
     9  -.251052453492112032869491235353891D-18,
     *  0.118067019426474278738319173650695D-18,
     1  0.182746800335155441155764424855902D-18,
     2  -.642344360403561579751828672744009D-19/
c
        save
c
        n=42
        done=1
        gam=0
c
        t=(x-1-done/2)*2 
c
        tt=1
        do 1200 i=1,n
c
        gam=gam+coefs(i)*tt
        tt=tt*t
 1200 continue
c
        return
        end
c
c 
c 
c 
c 
        subroutine gammanew_eval2(x,gam) 
        implicit real *8 (a-h,o-z)
        dimension coefs(22)
c
        data coefs/
     1  0.886226925452758027D+00,0.161691987244425025D-01,
     2  0.103703363422071934D+00,-.134118505705954070D-01,
     3  0.904033494042846197D-02,-.242259538441698247D-02,
     4  0.915785994891075934D-03,-.296890120771614378D-03,
     5  0.100928168758020265D-03,-.336375903860729251D-04,
     6  0.112523678860543027D-04,-.375495650035445620D-05,
     7  0.125310464830949639D-05,-.417909902825316389D-06,
     8  0.138824572934630909D-06,-.462835109086989667D-07,
     9  0.160743202007499598D-07,-.535845567983428334D-08,
     *  0.129318662203785113D-08,-.431075842242065260D-09,
     1  0.356478271052638709D-09,-.118826785520470270D-09/
        save
c
        n=22
        done=1
        gam=0
c
        t=(x-1-done/2)*2 
c
        tt=1
        do 1200 i=1,n
c
        gam=gam+coefs(i)*tt
        tt=tt*t
 1200 continue
c
        return
        end

