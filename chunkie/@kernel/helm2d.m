function obj = helm2d(type, zk, coefs)
%KERNEL.HELM2D   Construct the Helmholtz kernel.
%   KERNEL.HELM2D('s', ZK) or KERNEL.HELM2D('single', ZK) constructs the
%   single-layer Helmholtz kernel with wavenumber ZK.
%
%   KERNEL.HELM2D('d', ZK) or KERNEL.HELM2D('double', ZK) constructs the
%   double-layer Helmholtz kernel with wavenumber ZK.
%
%   KERNEL.HELM2D('sp', ZK) or KERNEL.HELM2D('sprime', ZK) constructs the
%   derivative of the single-layer Helmholtz kernel with wavenumber ZK.
%
%   KERNEL.HELM2D('c', ZK, COEFS) or KERNEL.HELM2D('combined', ZK, COEFS)
%   constructs the combined-layer Helmholtz kernel with wavenumber ZK and
%   parameter ETA, i.e., COEFS(1)*KERNEL.HELM2D('d', ZK) + 
%   COEFS(2)*KERNEL.HELM2D('s', ZK).
%
%   KERNEL.HELM2D('cp', ZK, COEFS) or KERNEL.HELM2D('combined', ZK, COEFS)
%   constructs the derivative of the combined-layer Helmholtz kernel with 
%   wavenumber ZK and  parameter ETA, i.e., COEFS(1)*KERNEL.HELM2D('dp', ZK) + 
%   COEFS(2)*KERNEL.HELM2D('sp', ZK).
%
%   KERNEL.HELM2D('all', ZK, COEFS) or KERNEL.HELM2D('tsys', ZK, COEFS) 
%   constructs the (2x2) matrix of kernels [D, S; D' S'] scaled by coefs
%   where coefs is (2x2) matrix, i.e. D part is scaled as
%   COEFS(1,1)*KERNEL.HELM2D('d', zk), and so on.
%
%   KERNEL.HELM2D('trans_rep', ZK, COEFS) or KERNEL.HELM2D('trep', ZK, COEFS) 
%   constructs the transmission repretsentation, i.e. the (1x2) matrix of
%   kernels [D, S] scaled by coefs where coefs is (1x2) matrix, i.e. D part
%   is scaled as COEFS(1)*KERNEL.HELM2D('d', zk), and so on.
%
%   KERNEL.HELM2D('trans_rep_p', ZK, COEFS) or KERNEL.HELM2D('trep_p', ZK, COEFS) 
%   constructs the derivative of the transmission repretsentation, i.e. the
%   (1x2) matrix of kernels [D', S'] scaled by coefs where coefs is (1x2)
%   matrix, i.e. D part is scaled as COEFS(1)*KERNEL.HELM2D('dp', zk), and
%   so on.
%
%   KERNEL.HELM2D('c2trans', ZK, COEFS) or KERNEL.HELM2D('c2t', ZK, COEFS) 
%   evaluates the combined-layer Helmholtz kernel and its derivative, i.e.
%   the (2x1) matrix of kernels [C; C'] scaled by coefs where coefs is
%   (2x2) matrix, i.e. kernel returns 
%   [COEFS(1,1)*KERNEL.HELM2D('d',zk)+COEFS(1,2)*KERNEL.HELM2D('s', zk); 
%   COEFS(2,1)*KERNEL.HELM2D('dp',zk)+COEFS(2,2)*KERNEL.HELM2D('sp', zk)]
% See also CHNK.HELM2D.KERN.

% author: Dan Fortunato
  
if ( nargin < 1 )
    error('Missing Helmholtz kernel type.');
end

if ( nargin < 2 )
    error('Missing Helmholtz wavenumber.');
end

obj = kernel();
obj.name = 'helmholtz';
obj.params.zk = zk;
obj.opdims = [1 1];

switch lower(type)

    case {'s', 'single'}
        obj.type = 's';
        obj.eval = @(s,t) chnk.helm2d.kern(zk, s, t, 's');
        obj.fmm  = @(eps,s,t,sigma) chnk.helm2d.fmm(eps, zk, s, t, 's', sigma);
        obj.sing = 'log';
        obj.splitinfo = [];
        obj.splitinfo.type = {[0 0 0 0],[1 0 0 0]};
        obj.splitinfo.action = {'r','r'};
        obj.splitinfo.functions = @(s,t) helm2d_s_split(zk,s,t);

    case {'d', 'double'}
        obj.type = 'd';
        obj.eval = @(s,t) chnk.helm2d.kern(zk, s, t, 'd');
        obj.fmm  = @(eps,s,t,sigma) chnk.helm2d.fmm(eps, zk, s, t, 'd', sigma);
        obj.sing = 'log';
        obj.splitinfo = [];
        obj.splitinfo.type = {[0 0 0 0],[1 0 0 0],[0 0 -1 0]};
        obj.splitinfo.action = {'r','r','r'};
        obj.splitinfo.functions = @(s,t) helm2d_d_split(zk,s,t);

    case {'sp', 'sprime'}
        obj.type = 'sp';
        obj.eval = @(s,t) chnk.helm2d.kern(zk, s, t, 'sprime');
        obj.fmm  = @(eps,s,t,sigma) chnk.helm2d.fmm(eps, zk, s, t, 'sprime', sigma);
        obj.sing = 'log';

    case {'dp', 'dprime'}
        obj.type = 'dp';
        obj.eval = @(s,t) chnk.helm2d.kern(zk, s, t, 'dprime');
        obj.fmm  = @(eps,s,t,sigma) chnk.helm2d.fmm(eps, zk, s, t, 'dprime', sigma);
        obj.sing = 'hs';

    case {'c', 'combined'}
        if ( nargin < 3 )
            warning('Missing combined layer coefficients. Defaulting to [1,1i].');
            coefs = [1,1i];
        end
        obj.type = 'c';
        obj.params.coefs = coefs;
        obj.eval = @(s,t) chnk.helm2d.kern(zk, s, t, 'c', coefs);
        obj.fmm  = @(eps,s,t,sigma) chnk.helm2d.fmm(eps, zk, s, t, 'c', sigma, coefs);
        obj.sing = 'log';
        obj.splitinfo = [];
        obj.splitinfo.type = {[0 0 0 0],[1 0 0 0],[0 0 -1 0]};
        obj.splitinfo.action = {'r','r','r'};
        obj.splitinfo.functions = @(s,t) helm2d_c_split(zk,s,t,coefs);

    case {'cp', 'cprime'}
        if ( nargin < 3 )
            warning('Missing combined layer coefficients. Defaulting to [1,1i].');
            coefs = [1,1i];
        end
        obj.type = 'cprime';
        obj.params.coefs = coefs;
        obj.eval = @(s,t) chnk.helm2d.kern(zk, s, t, 'cprime', coefs);
        obj.fmm  = @(eps,s,t,sigma) chnk.helm2d.fmm(eps, zk, s, t, 'cprime', sigma, coefs);
        obj.sing = 'hs';

    case {'all', 'trans_sys', 'tsys'}
        if ( nargin < 3 )
            warning('Missing transmission coefficients. Defaulting to [1,1;1,1].');
            coefs = ones(2,2);
        end
        obj.type = 'all';
        obj.opdims = [2,2];
        obj.eval = @(s,t) chnk.helm2d.kern(zk, s, t, 'all', coefs);
        obj.fmm  = []; 
        obj.sing = 'hs'; 

    case {'trans_rep','trep'} 
        if ( nargin < 3 )
            warning('Missing transmission coefficients. Defaulting to [1;1].');
            coefs = ones(2,1);
        end
        obj.type = 'trep';
        obj.opdims = [1,2];
        obj.eval = @(s,t) chnk.helm2d.kern(zk, s, t, 'trans_rep', coefs);
        obj.fmm  = []; 
        obj.sing = 'pv'; 
    
    case {'trans_rep_prime','trep_p', 'trans_rep_p'}
        if ( nargin < 3 )
            warning('Missing transmission coefficients. Defaulting to [1;1].');
            coefs = ones(2,1);
        end
        obj.type = 'trep_p';
        obj.opdims = [1,2];
        obj.eval = @(s,t) chnk.helm2d.kern(zk, s, t, 'trans_rep_prime', coefs);
        obj.fmm  = []; 
        obj.sing = 'hs'; 

    case {'c2trans', 'c2t', 'c2tr'}
        if ( nargin < 3 )
            warning('Missing combined layer coefficients. Defaulting to [1,1i;1,1i].');
            coefs = [1,1i;1,1i];
        end
        obj.type = 'c2tr';
        obj.opdims = [2,1];
        obj.eval = @(s,t) chnk.helm2d.kern(zk, s, t, 'c2trans', coefs);
        obj.fmm  = []; 
        obj.sing = 'hs'; 
       
    otherwise
        error('Unknown Helmholtz kernel type ''%s''.', type);

end

icheck = exist(['fmm2d.' mexext], 'file');
if icheck ~=3
    obj.fmm = [];
end


end

function f = helm2d_s_split(zk,s,t)
seval = chnk.helm2d.kern(zk, s, t, 's');
dist = (s.r(1,:)+1i*s.r(2,:))-(t.r(1,:)'+1i*t.r(2,:)');
logeval = log(abs(dist));
f = cell(2,1);
f{1} = seval + 1/(2*pi)*4*logeval.*imag(seval);
f{2} = 4*imag(seval);
end

function f = helm2d_d_split(zk,s,t)
deval = chnk.helm2d.kern(zk, s, t, 'd');
dist = (s.r(1,:)+1i*s.r(2,:))-(t.r(1,:)'+1i*t.r(2,:)');
logeval = log(abs(dist));
cauchyeval = (s.n(1,:)+1i*s.n(2,:))./(dist);
f = cell(3, 1);
f{1} = deval + 1/(2*pi)*4*logeval.*imag(deval) + 1/(2*pi)*real(cauchyeval);
f{2} = 4*imag(deval);
f{3} = ones(size(t.r,2),size(s.r,2));
end

function f = helm2d_c_split(zk,s,t,coefs)
seval = chnk.helm2d.kern(zk, s, t, 's');
deval = chnk.helm2d.kern(zk, s, t, 'd');
dist = (s.r(1,:)+1i*s.r(2,:))-(t.r(1,:)'+1i*t.r(2,:)');
logeval = log(abs(dist));
cauchyeval = (s.n(1,:)+1i*s.n(2,:))./(dist);
f = cell(3, 1);
f{1} = coefs(1) * (deval + 2/pi*logeval.*imag(deval) + 1/(2*pi)*real(cauchyeval)) + ...
       coefs(2) * (seval + 2/pi*logeval.*imag(seval));
f{2} = coefs(1) * (4*imag(deval)) + coefs(2) * (4*imag(seval));
f{3} = coefs(1)*ones(size(t.r,2),size(s.r,2));
end

