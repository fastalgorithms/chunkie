function varargout = pquadwts(r,d,n,d2,wts,j,rt,t,w,opts,intp_ab,intp,types,ifup)
%CHNK.pquadwts product integration for interaction of kernel on chunk
% at targets
%
% WARNING: this routine is not designed to be user-callable and assumes
%   a lot of precomputed values as input
%
% Syntax: [varargout] = pquadwts(r,d,n,d2,wts,j,rt,t,w,opts,intp_ab,intp,types)
%
% Input:
%   r - chnkr nodes
%   d - chnkr derivatives at nodes
%   n - chnkr normals at nodes
%   d2 - chnkr 2nd derivatives at nodes
%   wts - chnkr integration weights for smooth functions
%   j - chunk of interest
%   rt - position of target points. if any are not used by kernel (or not
%        well defined, e.g. when not on curve), a dummy array of the
%        appropriate size should be supplied
%   t - (Legendre) integration nodes
%   w - (Legendre) integration weights
%   opts - opts.side
%   intp_ab - panel endpoints interpolation matrix
%   types - specified singularity types
%       [0 0  0 0] - smooth quadr
%       [1 0  0 0] - log singularity
%       [0 0 -1 0] - cauchy singularity
%       [0 0 -2 0] - hypersingular
%       [0 0 -3 0] - supersingular
%       each [a, b, c, d]
%       corresponds to (log(z-w))^a (log(zc-wc))^b (z-w)^c (zc-wc)^d
%
% Output
%   varargout - integration matrices for specified singularity types

if nargin < 14, ifup = false; end

% Helsing-Ojala (interior/exterior?)
xlohi = intp_ab*(r(1,:,j)'+1i*r(2,:,j)');         % panel end points
r_i = intp*(r(1,:,j)'+1i*r(2,:,j)');              % upsampled panel
d_i = (intp*(d(1,:,j)'+1i*d(2,:,j)'));            % r'
d2_i = (intp*(d2(1,:,j)'+1i*d2(2,:,j)'));         % r''
wts_i = wts(:,j)';                                %
sp = abs(d_i); tang = d_i./sp;                    % speed, tangent
n_i = -1i*tang;                                   % normal
cur = -real(conj(d2_i).*n_i)./sp.^2;              % curvature
wxp_i = w.*d_i;                                   % complex speed weights (Helsing's wzp)

% determine number of matrices needed
nout = 0;
for j = 1:length(types)
    type0 = types{j};
    if (all(type0 == [0 0 0 0]))
        nout = max(nout,0);
    elseif (all(type0 == [1 0 0 0]))
        nout = max(nout,1);
    elseif (all(type0 == [0 0 -1 0]))
        nout = max(nout,2);
    elseif (all(type0 == [0 0 -2 0]))
        nout = max(nout,3);
    elseif (all(type0 == [0 0 -3 0]))
        nout = max(nout,4);
    else
        error("Split panel quad type [%s] not available",...
            join(string(type0)," "));
    end
end

Ac = cell(4,1);
if nout > 0
    [Ac{1:nout}] = SDspecialquad(struct('x',rt(1,:)' + 1i*rt(2,:)'), ...
        struct('x',r_i,'nx',n_i,'wxp',wxp_i),xlohi(1),xlohi(2),opts.side);
end

if ifup
    wts_i = (w.*sp)';
else
    for j = 1:nout
        Ac{j} = Ac{j}*intp;
    end
end

varargout = cell(size(types));
for j = 1:length(types)
    type0 = types{j};
    if (all(type0 == [0 0 0 0]))
        varargout{j} = ones(size(rt,2),numel(wts_i)).*wts_i;
    elseif (all(type0 == [1 0 0 0]))
        varargout{j} = Ac{1};
    elseif (all(type0 == [0 0 -1 0]))
        varargout{j} = Ac{2};
    elseif (all(type0 == [0 0 -2 0]))
        varargout{j} = Ac{3};
    elseif (all(type0 == [0 0 -3 0]))
        varargout{j} = Ac{4};
    else
        error("Split panel quad type %3d not available", type0);
    end
end

end

function [As, A, A1, A2, A3, A4] = SDspecialquad(t,s,a,b,side)
% S+D together... (with Sn) should be the same as requesting S or D
% more expensive if request Dn...
% https://github.com/ahbarnett/BIE2D/blob/master/panels/LapDLP_closepanel.m

zsc = (b-a)/2; zmid = (b+a)/2; % rescaling factor and midpoint of src segment
y = (s.x-zmid)/zsc; x = (t.x-zmid)/zsc;  % transformed src nodes, targ pts
N = numel(x);                            % # of targets
p = numel(s.x);                          % assume panel order is # nodes
if N*p==0
    As=0; A=0; A1=0; A2=0;
    return
end
c = (1-(-1).^(1:p))./(1:p);              % Helsing c_k, k = 1..p.
V = ones(p,p); for k=2:p, V(:,k) = V(:,k-1).*y; end  % Vandermonde mat @ nodes
P = zeros(p+1,N);      % Build P, Helsing's p_k vectorized on all targs...
d = 1.1; inr = abs(x)<=d; ifr = abs(x)>d;      % near & far treat separately

gam = exp(1i*pi/4);  % smaller makes cut closer to panel. barnett 4/17/18
if side == 'e', gam = conj(gam); end   % note gam is a phase, rots branch cut
P(1,:) = log(gam) + log((1-x)./(gam*(-1-x)));  % init p_1 for all targs int

% upwards recurrence for near targets, faster + more acc than quadr...
% (note rotation of cut in log to -Im; so cut in x space is lower unit circle)
if N>1 || (N==1 && inr==1) % Criterion added by Bowei Wu 03/05/15 to ensure inr not empty
    for k=1:p
        P(k+1,inr) = x(inr).'.*P(k,inr) + c(k);
    end  % recursion for p_k
end
% fine quadr (no recurrence) for far targets (too inaccurate cf downwards)...
Nf = numel(find(ifr)); wxp = s.wxp/zsc; % rescaled complex speed weights
if Nf>0 % Criterion added by Bowei Wu 03/05/15 to ensure ifr not empty
    P(end,ifr) = sum(((wxp.*(V(:,end).*y(:)))*ones(1,Nf))./bsxfun(@minus,y,x(ifr).'));  % int y^p/(y-x)
    for ii = p:-1:2
        P(ii,ifr) = (P(ii+1,ifr)-c(ii))./x(ifr).';
    end
end

Q = zeros(p,N); % compute q_k from p_k via Helsing 2009 eqn (18)... (p even!)
% Note a rot ang appears here too...  4/17/18
%gam = exp(1i*pi/4); % 1i;  % moves a branch arc as in p_1
%if side == 'e', gam = conj(gam); end   % note gam is a phase, rots branch cut
Q(1:2:end,:) = P(2:2:end,:) - repmat(log((1-x.').*(-1-x.')),[ceil(p/2) 1]); % guessed!
% (-1)^k, k odd, note each log has branch cut in semicircle from -1 to 1
% Q(2:2:end,:) = P(3:2:end,:) - repmat(log((1-x.')./((-1-x.'))),[floor(p/2) 1]);  % same cut as for p_1
Q(2:2:end,:) = P(3:2:end,:) - repmat(log(gam) + log((1-x.')./(gam*(-1-x.'))),[floor(p/2) 1]);  % same cut as for p_1
Q = Q.*repmat(1./(1:p)',[1 N]); % k even
As = real((V.'\Q).'.*repmat((1i*s.nx)',[N 1])*zsc)/(2*pi*abs(zsc));
As = As*abs(zsc) - log(abs(zsc))/(2*pi)*repmat(abs(s.wxp)',[N 1]); % unscale, yuk
warning('off','MATLAB:nearlySingularMatrix');
if nargout > 1
    % A = real((V.'\P).'*(1i/(2*pi)));         % solve for special quadr weights
    A = ((V.'\P(1:p,:)).'*(1i/(2*pi)));         % do not take real for the eval of Stokes DLP non-laplace term. Bowei 10/19/14
    %A = (P.'*inv(V))*(1i/(2*pi));   % equiv in exact arith, but not bkw stable.
end
if nargout > 2
    R =  -(kron(ones(p,1),1./(1-x.')) + kron((-1).^(0:p-1).',1./(1+x.'))) +...
        repmat((0:p-1)',[1 N]).*[zeros(1,N); P(1:p-1,:)];  % hypersingular kernel weights of Helsing 2009 eqn (14)
    Az = (V.'\R).'*(1i/(2*pi*zsc));  % solve for targ complex-deriv mat & rescale
    A1 = Az;
    if nargout > 3
        S = -(kron(ones(p,1),1./(1-x.').^2) - kron((-1).^(0:p-1).',1./(1+x.').^2))/2 +...
            repmat((0:p-1)',[1 N]).*[zeros(1,N); R(1:p-1,:)]/2; % supersingular kernel weights
        Azz = (V.'\S).'*(1i/(2*pi*zsc^2));
        if nargout > 4
            A1 = real(Az); A2 = -imag(Az);  % note sign for y-deriv from C-deriv
            A3 = real(Azz); A4 = -imag(Azz);
        else
            A1 = Az; A2 = Azz;
        end
    end
end

end
