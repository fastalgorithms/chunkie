function obj = stok2d(type, mu, coefs)
%KERNEL.STOK2D   Construct the Stokes kernel.
%   KERNEL.STOK2D('svel', MU) or KERNEL.ELAST2D('svelocity', MU)
%   constructs the single-layer Stokes kernel for velocity with viscosity
%   MU. KERNEL.STOK2D('s', MU) and KERNEL.STOK2D('single', MU) are
%   equivalent.
%
%   KERNEL.STOK2D('spres', MU) or KERNEL.ELAST2D('spressure', MU)
%   constructs the single-layer Stokes kernel for pressure with viscosity
%   MU.
%
%   KERNEL.STOK2D('strac', MU) or KERNEL.ELAST2D('straction', MU)
%   constructs the single-layer Stokes kernel for traction with viscosity
%   MU.
%
%   KERNEL.STOK2D('sgrad', MU) or KERNEL.ELAST2D('sgradient', MU)
%   constructs the single-layer Stokes kernel for velocity gradient with 
%   viscosity MU.
%
%   KERNEL.STOK2D('dvel', MU) or KERNEL.ELAST2D('dvelocity', MU)
%   constructs the double-layer Stokes kernel for velocity with viscosity
%   MU. KERNEL.STOK2D('d', MU) and KERNEL.STOK2D('double', MU) are
%   equivalent.
%
%   KERNEL.STOK2D('dpres', MU) or KERNEL.ELAST2D('dpressure', MU)
%   constructs the double-layer Stokes kernel for pressure with viscosity
%   MU.
%
%   KERNEL.STOK2D('dtrac', MU) or KERNEL.ELAST2D('dtraction', MU)
%   constructs the double-layer Stokes kernel for traction with viscosity
%   MU.
%
%   KERNEL.STOK2D('dgrad', MU) or KERNEL.ELAST2D('dgradient', MU)
%   constructs the double-layer Stokes kernel for velocity gradient with 
%   viscosity MU.
%  
%   KERNEL.STOK2D('cvel', MU) or KERNEL.ELAST2D('cvelocity', MU)
%   constructs the combined-layer Stokes kernel for velocity with viscosity
%   MU. KERNEL.STOK2D('c', MU) and KERNEL.STOK2D('comb', MU) are
%   equivalent.
%
%   KERNEL.STOK2D('cpres', MU) or KERNEL.ELAST2D('cpressure', MU)
%   constructs the combined-layer Stokes kernel for pressure with viscosity
%   MU.
%
%   KERNEL.STOK2D('ctrac', MU) or KERNEL.ELAST2D('ctraction', MU)
%   constructs the combined-layer Stokes kernel for traction with viscosity
%   MU.
%
%   KERNEL.STOK2D('cgrad', MU) or KERNEL.ELAST2D('cgradient', MU)
%   constructs the combined-layer Stokes kernel for velocity gradient with 
%   viscosity MU.
%
% See also CHNK.STOK2D.KERN.

% author: Dan Fortunato  

if ( nargin < 1 )
    error('Missing Stokes kernel type.');
end

if ( nargin < 2 )
    error('Missing Stokes viscosity mu.');
end

obj = kernel();
obj.name = 'stokes';
obj.params.mu = mu;

% TODO: Add FMMs.
% TODO: Add singularities for Stokes.

switch lower(type)

    case {'svel', 'svelocity', 's', 'single'}
        obj.type = 'svel';
        obj.eval = @(s,t) chnk.stok2d.kern(mu, s, t, 'svel');
        obj.fmm = @(eps, s, t, sigma) chnk.stok2d.fmm(eps, mu, s, t, 'svel', sigma);
        obj.opdims = [2, 2];
        obj.sing = 'log';
        obj.splitinfo = [];
        obj.splitinfo.type = {[1 0 0 0],[0 0 -1 0],[0 0 -1 0]};
        obj.splitinfo.action = {'r','r','i'};
        obj.splitinfo.functions = @(s,t) stok2d_s_split(mu, s, t);

    case {'spres', 'spressure'}
        obj.type = 'spres';
        obj.eval = @(s,t) chnk.stok2d.kern(mu, s, t, 'spres');
        obj.fmm = @(eps, s, t, sigma) chnk.stok2d.fmm(eps, mu, s, t, 'spres', sigma);
        obj.opdims = [1, 2];
	      obj.sing = 'pv';
        obj.splitinfo = [];
        obj.splitinfo.type = {[0 0 -1 0],[0 0 -1 0]};
        obj.splitinfo.action = {'r','i'};
        obj.splitinfo.functions = @(s,t) stok2d_spres_split(mu, s, t);

    case {'strac', 'straction'}
        obj.type = 'strac';
        obj.eval = @(s,t) chnk.stok2d.kern(mu, s, t, 'strac');
        obj.fmm = @(eps, s, t, sigma) chnk.stok2d.fmm(eps, mu, s, t, 'strac', sigma);
        obj.opdims = [2, 2];
        obj.sing = 'smooth';
        obj.splitinfo = [];
        obj.splitinfo.type = {[0 0 -1 0],[0 0 -1 0],[0 0 -2 0],[0 0 -2 0]};
        obj.splitinfo.action = {'r','i','r','i'};
        obj.splitinfo.functions = @(s,t) stok2d_strac_split(mu, s, t);

    case {'sgrad', 'sgradient'}
        obj.type = 'sgrad';
        obj.eval = @(s,t) chnk.stok2d.kern(mu, s, t, 'sgrad');
        obj.fmm = @(eps, s, t, sigma) chnk.stok2d.fmm(eps, mu, s, t, 'sgrad', sigma);
        obj.opdims = [4, 2];
        obj.sing = 'pv';    

    case {'dvel', 'dvelocity', 'd', 'double'}
        obj.type = 'dvel';
        obj.eval = @(s,t) chnk.stok2d.kern(mu, s, t, 'dvel');
        obj.fmm = @(eps, s, t, sigma) chnk.stok2d.fmm(eps, mu, s, t, 'dvel', sigma);
        obj.opdims = [2, 2];
        obj.sing = 'smooth'; 
        obj.splitinfo = [];
        obj.splitinfo.type = {[0 0 -1 0],[0 0 -1 0],[0 0 -2 0],[0 0 -2 0]};
        obj.splitinfo.action = {'r','i','r','i'};
        obj.splitinfo.functions = @(s,t) stok2d_d_split(mu, s, t);

    case {'dpres', 'dpressure'}
        obj.type = 'dpres';
        obj.eval = @(s,t) chnk.stok2d.kern(mu, s, t, 'dpres');
        obj.fmm = @(eps, s, t, sigma) chnk.stok2d.fmm(eps, mu, s, t, 'dpres', sigma);
        obj.opdims = [1, 2];
	      obj.sing = 'hs';
        obj.splitinfo = [];
        obj.splitinfo.type = {[0 0 -2 0],[0 0 -2 0]};
        obj.splitinfo.action = {'r','i'};
        obj.splitinfo.functions = @(s,t) stok2d_dpres_split(mu, s, t);

    case {'dtrac', 'dtraction'}
        obj.type = 'dtrac';
        obj.eval = @(s,t) chnk.stok2d.kern(mu, s, t, 'dtrac');
        obj.fmm = @(eps, s, t, sigma) chnk.stok2d.fmm(eps, mu, s, t, 'dtrac', sigma);
        obj.opdims = [2, 2];
        obj.sing = 'hs';
        obj.splitinfo = [];
        obj.splitinfo.type = {[0 0 -2 0],[0 0 -2 0],[0 0 -3 0],[0 0 -3 0]};
        obj.splitinfo.action = {'r','i','r','i'};
        obj.splitinfo.functions = @(s,t) stok2d_dtrac_split(mu, s, t);

    case {'dgrad', 'dgradient'}
        obj.type = 'dgrad';
        obj.eval = @(s,t) chnk.stok2d.kern(mu, s, t, 'dgrad');
        obj.fmm = @(eps, s, t, sigma) chnk.stok2d.fmm(eps, mu, s, t, 'dgrad', sigma);
        obj.opdims = [4, 2];
        obj.sing = 'hs';    


    case {'cvel', 'cvelocity', 'c', 'comb'}
        if ( nargin < 2 )
            warning('Missing combined layer parameter coefs. Defaulting to 1.');
            coefs = ones(2,1);
        end
        obj.type = 'cvel';
        obj.params.coefs = coefs;
        obj.eval = @(s,t) coefs(1)*chnk.stok2d.kern(mu, s, t, 'dvel') + ...
                          coefs(2)*chnk.stok2d.kern(mu, s, t, 'svel');
        obj.fmm = @(eps, s, t, sigma) chnk.stok2d.fmm(eps, mu, s, t, 'cvel', sigma, coefs);
        obj.opdims = [2, 2];
        obj.sing = 'log';

    case {'cpres', 'cpressure'}
        if ( nargin < 2 )
            warning('Missing combined layer parameter coefs. Defaulting to 1.');
            coefs = ones(2,1);
        end
        obj.type = 'cpres';
        obj.params.coefs = coefs;
        obj.eval = @(s,t) coefs(1)*chnk.stok2d.kern(mu, s, t, 'dpres') + ...
                          coefs(2)*chnk.stok2d.kern(mu, s, t, 'spres');
        obj.fmm = @(eps, s, t, sigma) chnk.stok2d.fmm(eps, mu, s, t, 'cpres', sigma, coefs);
        obj.opdims = [1, 2];
	    obj.sing = 'hs';

    case {'ctrac', 'ctraction'}
        if ( nargin < 2 )
            warning('Missing combined layer parameter coefs. Defaulting to 1.');
            coefs = ones(2,1);
        end
        obj.type = 'ctrac';
        obj.params.coefs = coefs;
        obj.eval = @(s,t) coefs(1)*chnk.stok2d.kern(mu, s, t, 'dtrac') + ...
                          coefs(2)*chnk.stok2d.kern(mu, s, t, 'strac');
        obj.fmm = @(eps, s, t, sigma) chnk.stok2d.fmm(eps, mu, s, t, 'ctrac', sigma, coefs);
        obj.opdims = [2, 2];
        obj.sing = 'hs';

    case {'cgrad', 'cgradient'}
        obj.type = 'cgrad';
       obj.eval = @(s,t) coefs(1)*chnk.stok2d.kern(mu, s, t, 'dgrad') + ...
                          coefs(2)*chnk.stok2d.kern(mu, s, t, 'sgrad');
        obj.fmm = @(eps, s, t, sigma) chnk.stok2d.fmm(eps, mu, s, t, 'cgrad', sigma, coefs);
        obj.opdims = [4, 2];
        obj.sing = 'hs';    

    otherwise
        error('Unknown Stokes kernel type ''%s''.', type);

end

icheck = exist(['fmm2d.' mexext], 'file');
if icheck ~=3
    obj.fmm = [];
end

end

function f = stok2d_s_split(mu,s,t)
dist = (s.r(1,:)+1i*s.r(2,:))-(t.r(1,:)'+1i*t.r(2,:)');
distr = real(dist);
disti = imag(dist);
f = cell(3, 1);
ntarg = numel(t.r(1,:));
nsrc = numel(s.r(1,:));
sn1mat = repmat(s.n(1,:),[ntarg 1]);
sn2mat = repmat(s.n(2,:),[ntarg 1]);
f{1} = zeros(2*ntarg,2*nsrc);
f{2} = zeros(2*ntarg,2*nsrc);
f{3} = zeros(2*ntarg,2*nsrc);
f{1}(1:2:end,1:2:end) = 1/(2*mu);
f{1}(2:2:end,2:2:end) = 1/(2*mu);
f{2}(1:2:end,1:2:end) = -sn1mat.*distr/(2*mu);
f{2}(1:2:end,2:2:end) = -sn2mat.*distr/(2*mu);
f{2}(2:2:end,1:2:end) = -sn1mat.*disti/(2*mu);
f{2}(2:2:end,2:2:end) = -sn2mat.*disti/(2*mu);
f{3}(1:2:end,1:2:end) = -sn2mat.*distr/(2*mu);
f{3}(1:2:end,2:2:end) =  sn1mat.*distr/(2*mu);
f{3}(2:2:end,1:2:end) = -sn2mat.*disti/(2*mu);
f{3}(2:2:end,2:2:end) =  sn1mat.*disti/(2*mu);
end

function f = stok2d_d_split(mu,s,t)
dist = (s.r(1,:)+1i*s.r(2,:))-(t.r(1,:)'+1i*t.r(2,:)');
distr = real(dist);
disti = imag(dist);
f = cell(4, 1);
ntarg = numel(t.r(1,:));
nsrc = numel(s.r(1,:));
snci = 1./(s.n(1,:)+1i*s.n(2,:));
sn11mat = repmat(real(snci).*s.n(1,:),[ntarg 1]);
sn12mat = repmat(real(snci).*s.n(2,:),[ntarg 1]);
sn21mat = repmat(imag(snci).*s.n(1,:),[ntarg 1]);
sn22mat = repmat(imag(snci).*s.n(2,:),[ntarg 1]);
f{1} = zeros(2*ntarg,2*nsrc);
f{2} = zeros(2*ntarg,2*nsrc);
f{3} = zeros(2*ntarg,2*nsrc);
f{4} = zeros(2*ntarg,2*nsrc);
f{1}(1:2:end,1:2:end) =  sn11mat;
f{1}(1:2:end,2:2:end) = -sn21mat;
f{1}(2:2:end,1:2:end) =  sn12mat;
f{1}(2:2:end,2:2:end) = -sn22mat;
f{2}(1:2:end,1:2:end) = -sn21mat;
f{2}(1:2:end,2:2:end) = -sn11mat;
f{2}(2:2:end,1:2:end) = -sn22mat;
f{2}(2:2:end,2:2:end) = -sn12mat;
f{3}(1:2:end,1:2:end) =  distr;
f{3}(1:2:end,2:2:end) =  disti;
f{4}(2:2:end,1:2:end) = -distr;
f{4}(2:2:end,2:2:end) = -disti;
end

function f = stok2d_strac_split(mu,s,t)
dist = (s.r(1,:)+1i*s.r(2,:))-(t.r(1,:)'+1i*t.r(2,:)');
distr = real(dist);
disti = imag(dist);
f = cell(4, 1);
ntarg = numel(t.r(1,:));
nsrc = numel(s.r(1,:));
snci = 1./(s.n(1,:)+1i*s.n(2,:));
tndotsn = ((t.n(1,:)'+1i*t.n(2,:)').*((-1i*s.n(1,:)-s.n(2,:))/1i));
sndotdistdottn = ((ones(ntarg,1)*snci).*real((-distr+1i*disti).*(t.n(1,:)'+1i*t.n(2,:)')));
f{1} = zeros(2*ntarg,2*nsrc);
f{2} = zeros(2*ntarg,2*nsrc);
f{3} = zeros(2*ntarg,2*nsrc);
f{4} = zeros(2*ntarg,2*nsrc);
f{1}(1:2:end,1:2:end) =  real(tndotsn);
f{1}(2:2:end,2:2:end) =  real(tndotsn);
f{2}(1:2:end,1:2:end) = -imag(tndotsn);
f{2}(2:2:end,2:2:end) = -imag(tndotsn);
f{3}(1:2:end,1:2:end) =  real(sndotdistdottn);
f{3}(1:2:end,2:2:end) = -imag(sndotdistdottn);
f{3}(2:2:end,1:2:end) = -imag(sndotdistdottn);
f{3}(2:2:end,2:2:end) = -real(sndotdistdottn);
f{4}(1:2:end,1:2:end) = -imag(sndotdistdottn);
f{4}(1:2:end,2:2:end) = -real(sndotdistdottn);
f{4}(2:2:end,1:2:end) = -real(sndotdistdottn);
f{4}(2:2:end,2:2:end) =  imag(sndotdistdottn);
end

function f = stok2d_dtrac_split(mu,s,t)
dist = (s.r(1,:)+1i*s.r(2,:))-(t.r(1,:)'+1i*t.r(2,:)');
ntarg = numel(t.r(1,:));
nsrc = numel(s.r(1,:));
tn1mat = repmat(t.n(1,:)',[1 nsrc]);
tn2mat = repmat(t.n(2,:)',[1 nsrc]);
conjsncmat = repmat((s.n(1,:)-1i*s.n(2,:))./(s.n(1,:)+1i*s.n(2,:)),[ntarg 1]);
f = cell(4, 1);
f{1} = zeros(2*ntarg,2*nsrc);
f{2} = zeros(2*ntarg,2*nsrc);
f{3} = zeros(2*ntarg,2*nsrc);
f{4} = zeros(2*ntarg,2*nsrc);
f{1}(1:2:end,1:2:end) =    tn1mat + real(conjsncmat).*tn1mat - imag(conjsncmat).*tn2mat;
f{1}(1:2:end,2:2:end) =   -tn2mat - real(conjsncmat).*tn2mat - imag(conjsncmat).*tn1mat;
f{1}(2:2:end,1:2:end) =   -tn2mat - real(conjsncmat).*tn2mat - imag(conjsncmat).*tn1mat;
f{1}(2:2:end,2:2:end) =  3*tn1mat - real(conjsncmat).*tn1mat + imag(conjsncmat).*tn2mat;
f{1} = mu*f{1};
f{2}(1:2:end,1:2:end) = -3*tn2mat - real(conjsncmat).*tn2mat - imag(conjsncmat).*tn1mat;
f{2}(1:2:end,2:2:end) =    tn1mat - real(conjsncmat).*tn1mat + imag(conjsncmat).*tn2mat;
f{2}(2:2:end,1:2:end) =    tn1mat - real(conjsncmat).*tn1mat + imag(conjsncmat).*tn2mat;
f{2}(2:2:end,2:2:end) =   -tn2mat + real(conjsncmat).*tn2mat + imag(conjsncmat).*tn1mat;
f{2} = mu*f{2};
f{3}(1:2:end,1:2:end) = -4*real(-dist.*(tn1mat - 1i*tn2mat));
f{3}(2:2:end,2:2:end) =  4*real(-dist.*(tn1mat - 1i*tn2mat));
f{3} = mu*f{3};
f{4}(1:2:end,2:2:end) =  4*real(-dist.*(tn1mat - 1i*tn2mat));
f{4}(2:2:end,1:2:end) =  4*real(-dist.*(tn1mat - 1i*tn2mat));
f{4} = mu*f{4};
end

function f = stok2d_spres_split(mu,s,t)
f = cell(2, 1);
ntarg = numel(t.r(1,:));
nsrc = numel(s.r(1,:));
sn1mat = repmat(s.n(1,:),[ntarg 1]);
sn2mat = repmat(s.n(2,:),[ntarg 1]);
f{1} = zeros(ntarg,2*nsrc);
f{2} = zeros(ntarg,2*nsrc);
f{1}(:,1:2:end) =  sn1mat;
f{1}(:,2:2:end) =  sn2mat;
f{2}(:,1:2:end) =  sn2mat;
f{2}(:,2:2:end) = -sn1mat;
end

function f = stok2d_dpres_split(mu,s,t)
f = cell(2, 1);
ntarg = numel(t.r(1,:));
nsrc = numel(s.r(1,:));
f{1} = zeros(ntarg,2*nsrc);
f{2} = zeros(ntarg,2*nsrc);
f{1}(:,1:2:end) = -2*mu;
f{2}(:,2:2:end) =  2*mu;
end
