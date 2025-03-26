function fints = chunkerkerneval_smooth(chnkr,kern,opdims,dens, ...
    targinfo,flag,opts)

if isa(kern,'kernel')
    kerneval = kern.eval;
else
    kerneval = kern;
end

flam = false;
accel = true;
forcefmm = false;

if nargin < 6
    flag = [];
end
if nargin < 7
    opts = [];
end
if isfield(opts,'flam'); flam = opts.flam; end
if isfield(opts,'accel'); accel = opts.accel; end
if isfield(opts,'forcefmm'); forcefmm = opts.forcefmm; end

k = chnkr.k;
nch = chnkr.nch;

assert(numel(dens) == opdims(2)*k*nch,'dens not of appropriate size')
dens = reshape(dens,opdims(2),k,nch);

w = chnkr.wstor;
[~,nt] = size(targinfo.r);

fints = zeros(opdims(1)*nt,1);

% assume smooth weights are good enough

% Sequence of checks, first see ifflam is set as it supercedes
% everything, if not flam, then check to see if the fmm
% exists and whether it should be used
% The number of sources set to 200 is currently a hack, 
% must be set based on opdims, accuracy, and kernel type
% considerations

imethod = 'direct';
if flam
    imethod = 'flam';
elseif isa(kern,'kernel') && ~isempty(kern.fmm)
    if forcefmm
        imethod = 'fmm';
    elseif accel
        if nt > 200 || chnkr.npt > 200
            imethod = 'fmm';
         end
     end
end

diamsrc = max(abs(chnkr.r(:)));
diamtarg = max(abs(targinfo.r(:)));
diam = max(diamsrc, diamtarg);

% flag targets that are within 1e-14 of sources
% these are automatically ignored by the fmm and 
% in chunkermat corrections
if ~(strcmpi(imethod,'fmm') && isempty(flag))
    flagslf = chnk.flagself(targinfo.r, chnkr.r, 1e-14*diam);
    if isempty(flagslf)
        selfzero = sparse(opdims(1)*size(targinfo.r(:,:),2), ...
            opdims(2)*chnkr.npt);
    else
        tmp = repmat((1:opdims(1))',opdims(2),1) + opdims(1)*(flagslf(1,:)-1);
        flagslftarg = tmp(:);
        tmp = repmat((1:opdims(2)),opdims(1),1);
        tmp = tmp(:) + opdims(2)*(flagslf(2,:)-1);
        flagslfsrc = tmp(:);
        
        selfzero = sparse(flagslftarg,flagslfsrc, 1e-300, ...
            opdims(1)*size(targinfo.r(:,:),2), opdims(2)*chnkr.npt);
    end
else
    selfzero = sparse(opdims(1)*size(targinfo.r(:,:),2), ...
        opdims(2)*chnkr.npt);
end

if strcmpi(imethod,'direct')
    % do dense version
    if isempty(flag)
        % nothing to ignore
        for i = 1:nch
            densvals = dens(:,:,i); densvals = densvals(:);
            dsdtdt = sqrt(sum(chnkr.d(:,:,i).^2,1));
            dsdtdt = dsdtdt(:).*w(:);
            dsdtdt = repmat( (dsdtdt(:)).',opdims(2),1);
            densvals = densvals.*(dsdtdt(:));
            srcinfo = []; srcinfo.r = chnkr.r(:,:,i); 
            srcinfo.n = chnkr.n(:,:,i);
            srcinfo.d = chnkr.d(:,:,i); srcinfo.d2 = chnkr.d2(:,:,i);
            kernmat = kerneval(srcinfo,targinfo);

            selfzeroch = selfzero(:, opdims(2)*k*(i-1) + (1:opdims(2)*k));
            [isp,jsp,~] = find(selfzeroch);
            linsp = isp + (jsp-1)*size(selfzeroch,1);
            kernmat(linsp) = 0;

            fints = fints + kernmat*densvals;
            % sum(fints)
        end
    else
        % ignore interactions in flag array
        for i = 1:nch
            densvals = dens(:,:,i); densvals = densvals(:);
            dsdtdt = sqrt(sum(chnkr.d(:,:,i).^2,1));
            dsdtdt = dsdtdt(:).*w(:);
            dsdtdt = repmat( (dsdtdt(:)).',opdims(2),1);
            densvals = densvals.*(dsdtdt(:));
            srcinfo = []; srcinfo.r = chnkr.r(:,:,i); 
            srcinfo.n = chnkr.n(:,:,i);
            srcinfo.d = chnkr.d(:,:,i); srcinfo.d2 = chnkr.d2(:,:,i);
            kernmat = kerneval(srcinfo,targinfo);

            rowkill = find(flag(:,i)); 
            rowkill = (opdims(1)*(rowkill(:)-1)).' + (1:opdims(1)).';
            kernmat(rowkill,:) = 0;

            selfzeroch = selfzero(:, opdims(2)*k*(i-1) + (1:opdims(2)*k));
            [isp,jsp,~] = find(selfzeroch);
            linsp = isp + (jsp-1)*size(selfzeroch,1);
            kernmat(linsp) = 0;

            fints = fints + kernmat*densvals;
        end
    end
else

    wts = chnkr.wts;
    wts = wts(:);
    
    if strcmpi(imethod,'flam')
        xflam1 = chnkr.r(:,:);
        xflam1 = repmat(xflam1,opdims(2),1);
        xflam1 = reshape(xflam1,chnkr.dim,numel(xflam1)/chnkr.dim);

        targinfo_flam = [];
        targinfo_flam.r = repelem(targinfo.r(:,:),1,opdims(1));
        if isfield(targinfo, 'd')
            targinfo_flam.d = repelem(targinfo.d(:,:),1,opdims(1));
        end
        
        if isfield(targinfo, 'd2')
            targinfo_flam.d2 = repelem(targinfo.d2(:,:),1,opdims(1));
        end
        
        if isfield(targinfo, 'n')
            targinfo_flam.n = repelem(targinfo.n(:,:),1,opdims(1));
        end

% TODO: Pull through data?

        matfun = @(i,j) chnk.flam.kernbyindexr(i, j, targinfo_flam, ...,
                           chnkr, kerneval, opdims, selfzero);
    

        width = max(abs(max(chnkr)-min(chnkr)))/3;
        tmax = max(targinfo.r(:,:),[],2); tmin = min(targinfo.r(:,:),[],2);
        wmax = max(abs(tmax-tmin));
        width = max(width,wmax/3);  
        npxy = chnk.flam.nproxy_square(kerneval,width);
        [pr,ptau,pw,pin] = chnk.flam.proxy_square_pts(npxy);

        verb = false; % TODO: make this an option to chunkerkerneval?
        optsnpxy = []; optsnpxy.rank_or_tol = opts.eps;
        pxyfun = @(lvl) proxyfunrbylevel(width,lvl,optsnpxy, ...
            chnkr,kerneval,opdims,verb && opts.proxybylevel);
        if ~opts.proxybylevel
            % if not using proxy-by-level, always use pxyfunr from level 1
            pxyfun = pxyfun(1);
        end

        optsifmm=[]; 
        optsifmm.Tmax=Inf; 
        optsifmm.proxybylevel = opts.proxybylevel;
        optsifmm.verb = verb;
        F = ifmm(matfun,targinfo_flam.r,xflam1,200,1e-14,pxyfun,optsifmm);
        fints = ifmm_mv(F,dens(:),matfun);
        fints = fints(:);
    else
        wts2 = repmat(wts(:).', opdims(2), 1);
        sigma = wts2(:).*dens(:);
        fints = kern.fmm(1e-14, chnkr, targinfo, sigma);
        fints = fints(:);
    end
    % delete interactions in flag array (possibly unstable approach)
    
    if ~isempty(flag)
        for i = 1:nch
            densvals = dens(:,:,i); densvals = densvals(:);
            dsdtdt = sqrt(sum(chnkr.d(:,:,i).^2,1));
            dsdtdt = dsdtdt(:).*w(:);
            dsdtdt = repmat( (dsdtdt(:)).',opdims(2),1);
            densvals = densvals.*(dsdtdt(:));
            srcinfo = []; srcinfo.r = chnkr.r(:,:,i); 
            srcinfo.n = chnkr.n(:,:,i);
            srcinfo.d = chnkr.d(:,:,i); srcinfo.d2 = chnkr.d2(:,:,i);

            delsmooth = find(flag(:,i)); 
            delsmoothrow = (opdims(1)*(delsmooth(:)-1)).' + (1:opdims(1)).';
            delsmoothrow = delsmoothrow(:);

            targinfo_use = [];
            targinfo_use.r = targinfo.r(:,delsmooth);

            if isfield(targinfo, 'd')
                targinfo_use.d = targinfo.d(:,delsmooth);
            end
        
            if isfield(targinfo, 'd2')
                targinfo_use.d2 = targinfo.d2(:,delsmooth);
            end
        
            if isfield(targinfo, 'n')
                targinfo_use.n = targinfo.n(:,delsmooth);
            end

            kernmat = kerneval(srcinfo,targinfo_use);

            selfzeroch = selfzero(:, opdims(2)*k*(i-1) + (1:opdims(2)*k));
            [isp,jsp,~] = find(selfzeroch);
            linsp = isp + (jsp-1)*length(i(:));
            kernmat(linsp) = 0;

            fints(delsmoothrow) = fints(delsmoothrow) - kernmat*densvals;
        end
    end    
end
% sum(fints)
end

function pxyfunrlvl = proxyfunrbylevel(width,lvl,optsnpxy, ...
    chnkr,kern,opdims,verb ...
    )
    npxy = chnk.flam.nproxy_square(kern,width/2^(lvl-1),optsnpxy);
    [pr,ptau,pw,pin] = chnk.flam.proxy_square_pts(npxy);

    if verb; fprintf('%3d | npxy = %i\n', lvl, npxy); end

    % return the FLAM-style pxyfunr corresponding to lvl
    pxyfunrlvl = @(rc,rx,cx,slf,nbr,l,ctr) chnk.flam.proxyfunr(rc,rx,slf,nbr,l, ...
        ctr,chnkr,kern,opdims,pr,ptau,pw,pin);
end