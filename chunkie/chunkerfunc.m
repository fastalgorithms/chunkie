function chnkr = chunkerfunc(fcurve,cparams,pref)
%CHUNKERFUNC create a chunker object corresponding to a parameterized curve
%
% Syntax: chnkr = chunkerfunc(fcurve,cparams,pref)
%
% Input: 
%   fcurve - function handle of the form
%               r = fcurve(t);
%            where r is a size [dim,size(t)] arrays describing
%            the position of a curve in dim dimensions parameterized by t.
%
%            optionally, the function can be of the form 
%               [r,d] = fcurve(t);  or [r,d,d2] = fcurve(t);
%            where d is the first derivative of r with respect to t and 
%            d2 is the second derivative. in some situations, this will
%            improve the convergence order and final precision.
%
% Optional input:
%	cparams - curve parameters structure (defaults)
%       cparams.ta = left end of t interval (0)
%       cparams.tb = right end of t interval (2*pi)
%       cparams.tsplits = set of initial break points for discretization in
%                 parameter space (should be in [ta,tb])
%       cparams.ifclosed = flag determining if the curve
%           is to be interpreted as a closed curve (true)
%       cparams.chsmall = max size of end intervals if
%           ifclosed == 0 (Inf)
%       cparams.nover = oversample resolved curve nover
%           times (0)
%       cparams.eps = tolerance to resolve coordinates and arclength 
%           density (1.0e-6)
%       cparams.lvlr = string, determines type of level
%          restriction to be enforced
%               lvlr = 'a' -> no chunk should have double the arc length 
%                               of its neighbor 
%               lvlr = 't' -> no chunk should have double the length in 
%                               parameter space of its neighbor 
%               lvlr = 'n' -> no enforcement of level restriction
%       cparams.lvlrfac = factor in level restriction, i.e. check if 
%               neighboring chunks differ in size by this factor (2.0)
%       cparams.maxchunklen - maximum length of any chunk (Inf)
%   pref - chunkerpref object or structure (defaults)
%       pref.nchmax - maximum number of chunks (10000)
%       pref.k - number of Legendre nodes on chunks (16)
%
% Output:
%   chnkr - a chunker object containing the discretization of the domain
%
% Examples:
%   chnkr = chunkerfunc(@(t) starfish(t)); % chunk up starfish w/ standard
%                                        % options
%   pref = []; pref.k = 30; 
%   cparams = []; cparams.eps = 1e-3;
%   chnkr = chunkerfunc(@(t) starfish(t),cparams,pref); % change up options
%   
% see also CHUNKERPOLY, CHUNKERPREF, CHUNKER

% author: Travis Askham (askhamwhat@gmail.com)
%

if nargin < 2
    cparams = [];
end
if nargin < 3
    pref = chunkerpref();
else
    pref = chunkerpref(pref);
end

chnkr = chunker(pref); % empty chunker

ta = 0.0; tb = 2*pi; ifclosed=true;
chsmall = Inf; nover = 0;
eps = 1.0e-6;
lvlr = 'a'; maxchunklen = Inf; lvlrfac = 2;
nout = 1;

if isfield(cparams,'ta')
    ta = cparams.ta;
end	 
if isfield(cparams,'tb')
    tb = cparams.tb;
end	 
if isfield(cparams,'ifclosed')
    ifclosed = cparams.ifclosed;
end	 
if isfield(cparams,'chsmall')
    chsmall = cparams.chsmall;
end	 
if isfield(cparams,'nover')
    nover = cparams.nover;
end	 
if isfield(cparams,'eps')
    eps = cparams.eps;
end	 
if isfield(cparams,'lvlr')
    lvlr = cparams.lvlr;
end
if isfield(cparams,'lvlrfac')
    lvlrfac = cparams.lvlrfac;
end
if isfield(cparams,'maxchunklen')
    maxchunklen = cparams.maxchunklen;
end

% discover number of outputs
try         
    [r,d,d2] = fcurve(t);
    nout = 3;
catch
    try 
        [r,d] = fcurve(ta);
        nout = 2;
    catch
        nout = 1;
    end
end

k = pref.k;
nchmax = pref.nchmax; 
 
dim = checkcurveparam(fcurve,ta,nout);
pref.dim = dim;
out = cell(3,1);

ifprocess = zeros(nchmax,1);

%       construct legendre nodes and weights, k and 2k of them, as well
%       as the interpolation/coefficients matrices

k2 = 2*k;
[xs,ws,us,vs] = lege.exps(k);
dermat = (vs*[lege.derpol(us); zeros(1,k)]).';
[xs2,ws2,u2,v2] = lege.exps(k2);   
dermat2 = (v2*[lege.derpol(u2); zeros(1,k2)]).';

xs2p = ((1:k2)-1)/(k2-1)*2-1;
[polvals,~] = lege.pols(xs2p,k-1);
interp_xs = reshape(polvals,[k,k2]).'*us; 

%       . . . start chunking

ab = zeros(2,nchmax);
adjs = zeros(2,nchmax);

if (isfield(cparams,'tsplits'))
    tsplits = cparams.tsplits;
    tsplits = [tsplits(:); ta; tb];
else
    tsplits = [ta;tb];
end
   
tsplits = sort(unique(tsplits),'ascend');
lab = length(tsplits);
if (lab-1 > nchmax)
    error(['CHUNKERFUNC: nchmax exceeded in chunkerfunc on initial splits.\n ',...
        'try increasing nchmax in preference struct']);
end
if (any(tsplits > tb) || any(tsplits < ta))
    error(['CHUNKERFUNC: tsplits outside of interval of definition.\n', ...
          'check definition of splits, ta and tb']);
end

ab(1,1:(lab-1)) = tsplits(1:end-1);
ab(2,1:(lab-1)) = tsplits(2:end);

nch=lab-1;
adjs(1,1:nch) = 0:(nch-1);
adjs(2,1:nch) = 2:(nch+1);

if ifclosed
    adjs(1,1)=nch;
    adjs(2,nch)=1;
else
    adjs(1,1)=-1;
    adjs(2,nch)=-1;
end
nchnew=nch;

maxiter_res=nchmax-nch;

rad_curr = 0;
for ijk = 1:maxiter_res

%       loop through all existing chunks, if resolved store, if not split
    xmin =  Inf;
    xmax = -Inf;
    ymin =  Inf;
    ymax = -Inf;
    
    ifdone=1;
    for ich=1:nchnew

        if (ifprocess(ich) ~= 1)
            ifprocess(ich)=1;

            a=ab(1,ich);
            b=ab(2,ich);
           
            ts = a + (b-a)*(xs2+1)/2.0;
            [out{1:nout}] = fcurve(ts);
            for j = nout+1:3
                out{j} = out{j-1}*dermat2*(2/(b-a));
            end
            r = out{1};
            d = out{2};
            d2 = out{3};
            zd = d(1,:)+1i*d(2,:);
            vd = abs(zd);
            zdd= d2(1,:)+1i*d2(2,:);
            dkappa = imag(zdd.*conj(zd))./abs(zd).^2;

            rlself = vd(:).'*ws2*(b-a)/2;
            
            cfs = u2*vd.';
            errs0 = sum(abs(cfs(1:k)).^2,1);
            errs = sum(abs(cfs(k+1:k2)).^2,1);
            err1 = sqrt(errs/errs0/k);
            
            resol_speed_test = err1>eps;
            if nout < 2
                resol_speed_test = err1>eps*k;
            end
            
            xmax = max(xmax,max(r(1,:)));
            ymax = max(ymax,max(r(2,:)));
            xmin = min(xmin,min(r(1,:)));
            ymin = min(ymin,min(r(2,:)));
            
            cfsx = u2*r(1,:).';
            cfsy = u2*r(2,:).';
            errx = sum(abs(cfsx(k+1:k2)).^2/k,1);
            erry = sum(abs(cfsy(k+1:k2)).^2/k,1);
            errx = sqrt(errx);
            erry = sqrt(erry);
            
            resol_curve_test = true;
            
            if (ijk >1)
                if (errx/rad_curr<eps && erry/rad_curr<eps)
                    resol_curve_test = false;
                end
            end    


           total_curve = (b-a)/2*sum(abs(dkappa).*ws2.');
           total_curve_test = total_curve >= (2*pi)/3;

           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            

     %       . . . mark as processed and resolved if less than eps

            if (resol_speed_test || resol_curve_test  || ... 
                   total_curve_test || rlself > maxchunklen || ...
                    and(or(adjs(1,ich) <= 0, adjs(2,ich) <= 0), ...
                rlself > chsmall))
              %       . . . if here, not resolved
              %       divide - first update the adjacency list
                if (nch +1 > nchmax)
                    error('CHUNKERFUNC: nchmax=%d exceeded. Unable to resolve curve.',nchmax)
                end

                ifprocess(ich)=0;
                ifdone=0;

                if ((nch == 1) && ifclosed)
                    adjs(1,nch)=2;
                    adjs(2,nch)=2;
                    adjs(1,nch+1)=1;
                    adjs(2,nch+1)=1;
                end

                if ((nch == 1) && (~ifclosed))
                    adjs(1,nch)=-1;
                    adjs(2,nch)=2;
                    adjs(1,nch+1)=1;
                    adjs(2,nch+1)=-1;
                end

                if (nch > 1)
                    iold2=adjs(2,ich);
                    adjs(2,ich)=nch+1;
                    if (iold2 > 0)
                        adjs(1,iold2)=nch+1;
                    end	
                    adjs(1,nch+1)=ich;
                    adjs(2,nch+1)=iold2;

                end
                    %       now update the endpoints in ab

                ab(1,ich)=a;
                ab(2,ich)=(a+b)/2;

                nch=nch+1;

                ab(1,nch)=(a+b)/2;
                ab(2,nch)=b;
            end
        end
    end
    if ((ifdone == 1) && (nchnew == nch))
        break;
    end
    nchnew=nch;
    
    rad_curr = max(xmax-xmin,ymax-ymin);
end


%       the curve should be resolved to precision eps now on
%       each interval ab(,i)
%       check the size of adjacent neighboring chunks - if off by a
%       factor of more than 2, split them as well. iterate until done.
   
if or(strcmpi(lvlr,'a'),strcmpi(lvlr,'t'))
    maxiter_adj=1000;
    for ijk = 1:maxiter_adj

        nchold=nch;
        ifdone=1;
        for i = 1:nchold
            i1=adjs(1,i);
            i2=adjs(2,i);

    %       calculate chunk lengths

            a=ab(1,i);
            b=ab(2,i);
            
            if strcmpi(lvlr,'a')
                rlself = chunklength(fcurve,a,b,xs,ws,nout,dermat);

                rl1=rlself;
                rl2=rlself;

                if (i1 > 0)
                    a1=ab(1,i1);
                    b1=ab(2,i1);
                    rl1 = chunklength(fcurve,a1,b1,xs,ws,nout,dermat);
                end
                if (i2 > 0)
                    a2=ab(1,i2);
                    b2=ab(2,i2);
                    rl2 = chunklength(fcurve,a2,b2,xs,ws,nout,dermat);
                end
            else
                
                rlself = b-a;
                rl1 = rlself;
                rl2 = rlself;
                if (i1 > 0)
                    rl1 = ab(2,i1)-ab(1,i1);
                end
                if (i2 > 0)
                    rl2 = ab(2,i2)-ab(1,i2);
                end
            end

    %       only check if self is larger than either of adjacent blocks,
    %       iterating a couple times will catch everything

            if (rlself > lvlrfac*rl1 || rlself > lvlrfac*rl2)

    %       split chunk i now, and recalculate nodes, ders, etc

                if (nch + 1 > nchmax)
                    error('CHUNKERFUNC: nchmax=%d exceeded during level restriction (curve resolved).',nchmax);
                end


                ifdone=0;
                a=ab(1,i);
                b=ab(2,i);
                ab2=(a+b)/2;

                i1=adjs(1,i);
                i2=adjs(2,i);
    %        
                adjs(1,i) = i1;
                adjs(2,i) = nch+1;

    %       . . . first update nch+1

                adjs(1,nch+1) = i;
                adjs(2,nch+1) = i2;

     %       . . . if there's an i2, update it

                if (i2 > 0)
                    adjs(1,i2) = nch+1;
                end

                nch=nch+1;
                ab(1,i)=a;
                ab(2,i)=ab2;

                ab(1,nch)=ab2;
                ab(2,nch)=b;
            end
        end

        if (ifdone == 1)
            break;
        end

    end
end

%       go ahead and oversample by nover, updating
%       the adjacency information adjs along the way


if (nover > 0) 
    for ijk = 1:nover

        nchold=nch;
        for i = 1:nchold
            a=ab(1,i);
            b=ab(2,i);

              % use dekker's method (hybrid secant/bisection)
           
            rl = chunklength(fcurve,a,b,xs,ws,nout,dermat);
            rlhalf=rl/2;
            
            fak = -rlhalf;
            fbk = rlhalf;
            
            thresh=eps*10*rad_curr;
            
            ak = a;
            bk = b;
            
            fbkm1 = fbk;
            bkm1 = bk;
            
            for iter = 1:52
                m = (ak+bk)/2;
                s = m;
                if fbkm1 ~= fbk
                    s = bk - (bk-bkm1)*fbk/(fbk-fbkm1);
                end
                
                if (s-m)*(s-bk) >= 0
                    s = m;
                end

                % store last iterate values
                bkm1 = bk;
                fbkm1 = fbk;

                bk = s;
                fbk = chunklength(fcurve,a,bk,xs,ws,nout,dermat)-rlhalf;
                
                % update bracket
                if fak*fbk >= 0
                    ak = bkm1;
                    fak = fbkm1;
                end
                
                % check if b is still the best guess
                if abs(fbk) > abs(fak)
                    tmp = bk;
                    bk = ak;
                    ak = tmp;
                    tmp = fbk;
                    fbk = fak;
                    fak = tmp;
                end                    
                
                if (abs(fbk) < thresh)
                    break
                    disp(iter)
                end
            end
	 
            ab2=bk;

            i1=adjs(1,i);
            i2=adjs(2,i);
            adjs(2,i)=nch+1;
            if (i2 > 0)
                adjs(1,i2)=nch+1;
            end

            if (nch + 1 > nchmax)
                error('CHUNKERFUNC: nchmax=%d exceeded while oversampling by nover=%d',nchmax,nover);
            end

            adjs(1,nch+1)=i;
            adjs(2,nch+1)=i2;
	 
            ab(1,i)=a;
            ab(2,i)=ab2;
	 
            nch=nch+1;

            ab(1,nch)=ab2;
            ab(2,nch)=b;
        end
    end
end

%       up to here, everything has been done in parameter space, [ta,tb]
%       . . . finally evaluate the k nodes on each chunk, along with 
%       derivatives and chunk lengths

chnkr = chnkr.addchunk(nch);


for i = 1:nch
    a=ab(1,i);
    b=ab(2,i);
    
    ts = a + (b-a)*(xs+1)/2;
    [out{1:nout}] = fcurve(ts);
    for j = nout+1:3
        out{j} = out{j-1}*dermat*(2/(b-a));
    end
    h = (b-a)/2;
    chnkr.rstor(:,:,i) = reshape(out{1},dim,k);
    chnkr.dstor(:,:,i) = reshape(out{2},dim,k)*h;
    chnkr.d2stor(:,:,i) = reshape(out{3},dim,k)*h*h;
end

chnkr.adjstor(:,1:nch) = adjs(:,1:nch);

% update normals
chnkr.nstor(:,:,1:nch) = normals(chnkr);

% update weights
chnkr.wtsstor(:,1:nch) = weights(chnkr);

end


function [len] = chunklength(fcurve,a,b,xs,ws,nout,dermat)
    
    out = cell(3,1);
    ts = a+(b-a)*(xs+1)/2;
    [out{1:nout}] = fcurve(ts);
    for j = nout+1:3
        out{j} = out{j-1}*dermat*(2/(b-a));
    end
    dsdt = sqrt(sum(abs(out{2}).^2,1));
    len = dot(dsdt,ws)*(b-a)/2;
 end

