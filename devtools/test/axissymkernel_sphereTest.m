iseed = 8675309;
rng(iseed);


zk = 40.1;

type = 'chnkr-star';

pref = [];
pref.k = 16;
ns = 10;
nt = 100;
ppw = 80;   % points per wavelength;
maxchunklen = pref.k/ppw/real(zk)*2*pi;
maxchunklen = 0.5;

[chnkr, ~, ~] = get_geometry(type, pref, ns, nt, maxchunklen);
wts = chnkr.wts; wts = wts(:);

fprintf('Done building geometry\n');

% Plot everything

figure(1)
clf
hold off
plot(chnkr)
axis equal

type = 'dprimediff';
K    = @(s,t) chnk.axissymhelm2d.kern(zk, s, t, [0,0], type);



% mode number
nmode = 20;

if strcmpi(type, 's')
    zfac = 1j*zk*sphericalbesselj(nmode,zk)*sphericalbesselh(nmode,zk);
elseif strcmpi(type, 'd')
    zfac = 1j*zk^2*sphericalbesselj_der(nmode,zk)*sphericalbesselh(nmode,zk) - 0.5;
elseif strcmpi(type, 'sprime')
    zfac = 1j*zk^2*sphericalbesselj(nmode,zk)*sphericalbesselh_der(nmode,zk) + 0.5;
elseif strcmpi(type, 'dprimediff')
    zfac = 1j*zk^3*sphericalbesselj_der(nmode,zk)*sphericalbesselh_der(nmode,zk);
    zfac = zfac - 1j*(1j*zk)^3*sphericalbesselj_der(nmode,1j*zk)*sphericalbesselh_der(nmode,1j*zk);
end



cos_t = chnkr.r(2,:)./sqrt(chnkr.r(2,:).^2 + chnkr.r(1,:).^2);
ubdry = legendre(nmode,chnkr.r(2,:));
ubdry = ubdry(1,:);
ubdry = ubdry(:);

% Form matrix
opts = [];
tic, A = chunkermat(chnkr, K, opts); toc

pot = A*ubdry;
potex = zfac*ubdry;

err1 = sqrt(sum(abs(pot - potex).^2.*wts(:)));
fprintf('error in layer pot on sphere=%d\n', err1);



function [chnkobj, sources, targets] = get_geometry(type, pref, ns, nt, maxchunklen)

if nargin == 0
    type = 'chnkr';
end

if nargin <= 1
    pref = [];
    pref.k = 16;
end

if nargin <= 2
    ns = 10;
end

if nargin <= 3
    nt = 10;
end


if nargin <= 4
    maxchunklen = 1.0;
end
    

if strcmpi(type, 'cgrph')
    
    
    nverts = 3; 
    verts = exp(-1i*pi/2 + 1i*pi*(0:(nverts-1))/(nverts-1));
    verts = [real(verts);imag(verts)];


    iind = 1:(nverts-1);
    jind = 1:(nverts-1);

    iind = [iind iind];
    jind = [jind jind + 1];
    jind(jind>nverts) = 1;
    svals = [-ones(1,nverts-1) ones(1,nverts-1)];
    edge2verts = sparse(iind, jind, svals, nverts-1, nverts);

    amp = 0.1;
    frq = 2;
    fchnks    = cell(1,size(edge2verts,1));
    for icurve = 1:size(edge2verts,1)
        fchnks{icurve} = @(t) sinearc(t, amp, frq);
    end
    cparams = [];
    cparams.nover = 2;
    cparams.maxchunklen = maxchunklen;

    chnkobj = chunkgraph(verts, edge2verts, fchnks, cparams, pref);
    chnkobj = balance(chnkobj);
    
    ts = 0.0+2*pi*rand(ns,1);
    sources = 3.0*[cos(ts)';sin(ts)'];
    
    ts = 0.0+2*pi*rand(nt,1);
    targets = 0.2*[cos(ts)'; sin(ts)'];


elseif strcmpi(type,'chnkr-star')
    cparams = [];
    cparams.eps = 1.0e-10;
    cparams.nover = 1;
    cparams.ifclosed = false;
    cparams.ta = -pi/2;
    cparams.tb = pi/2;
    cparams.maxchunklen = maxchunklen;
    narms = 0;
    amp = 0.0;
    chnkobj = chunkerfunc(@(t) starfish(t, narms, amp), cparams, pref); 
    chnkobj = sort(chnkobj);
    
    ts = -pi/2 + pi*rand(ns, 1);
    sources = starfish(ts, narms, amp);
    sources = 0.5*sources;

    ts = -pi/2 + pi*rand(nt, 1);
    targets = starfish(ts, narms, amp);
    targets = targets .* (1 + 0.5*repmat(rand(1, nt), 2, 1));   

elseif strcmpi(type,'chnkr-torus')
    cparams = [];
    cparams.eps = 1.0e-10;
    cparams.nover = 1;
    cparams.ifclosed = true;
    cparams.ta = 0;
    cparams.tb = 2*pi;
    cparams.maxchunklen = maxchunklen;
    narms = 0;
    amp = 0.25;
    ctr = [3 0];
    chnkobj = chunkerfunc(@(t) starfish(t, narms, amp, ctr), cparams, pref); 
    chnkobj = sort(chnkobj);
    
    ts = -pi/2 + pi*rand(ns, 1);
    sources = starfish(ts, narms, amp);
    sources = 0.5*sources;
    sources(1,:) = sources(1,:) + ctr(1);

    ts = -pi/2 + pi*rand(nt, 1);
    targets = starfish(ts, narms, amp);
    targets = targets .* (1 + 0.5*repmat(rand(1, nt), 2, 1));    
    targets(1,:) = targets(1,:) + ctr(1);
end



end


function [r,d,d2] = sinearc(t,amp,frq)
xs = t;
ys = amp*sin(frq*t);
xp = ones(size(t));
yp = amp*frq*cos(frq*t);
xpp = zeros(size(t));
ypp = -frq*frq*amp*sin(t);

r = [(xs(:)).'; (ys(:)).'];
d = [(xp(:)).'; (yp(:)).'];
d2 = [(xpp(:)).'; (ypp(:)).'];
end

function jn = sphericalbesselj(n,z)
    jn = sqrt(pi/2./z).*besselj(n+0.5,z);
end


function jnp = sphericalbesselj_der(n,z)
    jnp = 0.5*(sphericalbesselj(n-1,z) - (sphericalbesselj(n,z) + ...
                     z.*sphericalbesselj(n+1,z))./z);
end


function hn = sphericalbesselh(n,z)
    hn = sqrt(pi/2./z).*besselh(n+0.5,1,z);
end


function hnp = sphericalbesselh_der(n,z)
    hnp = 0.5*(sphericalbesselh(n-1,z) - (sphericalbesselh(n,z) + ...
                     z.*sphericalbesselh(n+1,z))./z);
end


