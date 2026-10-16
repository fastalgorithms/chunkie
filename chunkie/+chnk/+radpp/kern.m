function submat = kern(profun,dim,srcinfo,targinfo,type,varargin)
%CHNK.RADPP.KERN kernels of a radial Green's function in 2D or 3D
%
% Syntax: submat = chnk.radpp.kern(profun,dim,srcinfo,targinfo,type,varargin)
%
% The provided Green's function must be radial: G(x,y) = G(|x-y|) and
% return all the required derivatives in the usual format.
%
% S(x,y)  = G(x,y)
% D(x,y)  = G_{n_y}(x,y)
% S'(x,y) = G_{n_x}(x,y)
% D'(x,y) = G_{n_x n_y}(x,y)
%
% and so on.
%
% Input:
%   profun - function handle, profun(r) is the (numel(r),nprof) array of
%            on-axis profiles of G and its derivatives, see
%            CHNK.RADPP.PROFILES
%   dim - 2 or 3, the dimension of the derivative tensors, see
%            CHNK.RADPP.ROTATE
%   srcinfo - description of sources in ptinfo struct format, i.e.
%                ptinfo.r - positions (dim,:) array
%                ptinfo.d - first derivative in underlying
%                     parameterization (dim,:)
%                ptinfo.d2 - second derivative in underlying
%                     parameterization (dim,:)
%                ptinfo.n - normals (dim,:)
%                ptinfo.data - data (:,:), the first row is the
%                     derivative of the curvature for supported_eval
%   targinfo - description of targets in ptinfo struct format, only the
%                fields relevant to the kernel type are needed
%   type - string, determines kernel type, or a struct for a custom kernel
%                with fields fun, nderiv and opdims, see KERNEL.RADGREEN
%                type == 's', single layer kernel S
%                type == 'd', double layer kernel D
%                type == 'sp', normal derivative of single layer S'
%                type == 'dp', normal derivative of double layer D'
%                type == 'c', combined layer kernel coef(1) D + coef(2) S
%                type == 'cp', coef(1) D' + coef(2) S'
%                type == 'sg', 'dg', 'cg', gradients of S, D and
%                       coef(1) D + coef(2) S, a (dim*nt,ns) array
%                type == 'st', tangential derivative of S (dim = 2)
%                type == 'hilb', two times the adjoint of st (dim = 2)
%                type == 's2trans', [S; S']
%                type == 'd2trans', [D; D']
%                type == 'c2trans', returns the combined field and the
%                       normal derivative of the combined field
%                       [coef(1,1)*D + coef(1,2)*S; coef(2,1)*D' + coef(2,2)*S']
%                type == 'trans', returns all four layer potentials,
%                       [coef(1,1)*D coef(1,2)*S; coef(2,1)*D' coef(2,2)*S']
%                type == 'trans_rep', returns the potential corresponding
%                       to the transmission representation
%                       [coef(1)*D coef(2)*S]
%                type == 'trans_rep_prime', returns the normal derivative
%                       corresponding to the transmission representation
%                       [coef(1)*D' coef(2)*S']
%                type == 'trans_rep_grad', returns the gradient
%                       corresponding to the transmission representation
%                       [coef(1)*grad D coef(2)*grad S]
%
%         plate kernels (dim = 2), with n_x, tau_x and kappa_x the
%         normal, tangent and curvature of the targets and n_y, tau_y
%         and kappa_y those of the sources
%                type == 'clamped_bc', returns the clamped plate
%                       boundary conditions applied to a point source,
%                       [G; G_{n_x}]. requires target normal info.
%                type == 'clamped_eval', returns the kernels
%                       corresponding to the clamped plate representation,
%                       [G_{n_y n_y n_y}+3*G_{n_y tau_y tau_y},
%                        -G_{n_y n_y}+G_{tau_y tau_y}]
%                type == 'supported_bc', returns the supported plate
%                       boundary conditions applied to a point source,
%                       [G; G_{n_x n_x}+nu*G_{tau_x tau_x}]. requires
%                       target normal and tangent info.
%                type == 'supported_eval', returns the kernels
%                       corresponding to the supported plate
%                       representation. requires kappa' in the first
%                       data field of srcinfo.
%                type == 'free_bc', returns the free plate boundary
%                       conditions applied to a point source,
%                       [G_{n_x n_x}+nu*G_{tau_x tau_x};
%                        G_{n_x n_x n_x}+(2-nu)*G_{n_x tau_x tau_x}
%                         +(1-nu)*kappa_x*(G_{tau_x tau_x}-G_{n_x n_x})]
%                       requires target normals, tangents, and d2.
%                type == 'free_eval' returns the kernels corresponding to
%                       the free plate representation:
%                       [G_{n_y}, (1+nu)/2*G_{tau_y}, G]
%
%   varargin{1} - coef: length 2 array in the combined layer formulas,
%                 2x2 matrix for c2trans and trans, nu: scalar value for
%                 the Poisson's ratio of the plate in the supported and
%                 free plate kernels, otherwise does nothing
%
% Output:
%   submat - the evaluation of the selected kernel for the provided
%            sources and targets. the number of rows equals the number
%            of targets times opdims(1) and the number of columns equals
%            the number of sources times opdims(2), with the entries of
%            each block interleaved
%
% see also KERNEL.RADGREEN, KERNEL.RADPP, CHNK.RADPP.TYPEINFO

if isstruct(type)
    info = struct('type','custom','opdims',type.opdims, ...
        'nderiv',type.nderiv,'param','');
else
    info = chnk.radpp.typeinfo(type,dim);
end

par = [];
if ~isempty(varargin), par = varargin{1}; end
coef = [1 1];
nu = [];
switch info.param
    case 'coefs'
        if ~isempty(par), coef = par(:).'; end
    case 'coefs22'
        coef = ones(2);
        if ~isempty(par), coef = par; end
        if numel(coef) == 2, coef = repmat(coef(:).',2,1); end
    case 'nu'
        nu = par;
end

src = srcinfo.r(:,:);
targ = targinfo.r(:,:);

[~,ns] = size(src);
[~,nt] = size(targ);

dx = bsxfun(@minus,targ(1,:).',src(1,:));
dy = bsxfun(@minus,targ(2,:).',src(2,:));
dz = zeros(nt,ns);
if size(src,1) > 2 && size(targ,1) > 2
    dz = bsxfun(@minus,targ(3,:).',src(3,:));
end
r = sqrt(dx.^2 + dy.^2 + dz.^2);

% derivative tensors of G with respect to the target
comps = chnk.radpp.rotate(dx,dy,dz,profun(r(:)),info.nderiv,dim);
comps = reshape(comps,nt,ns,[]);

% derivative of G along the directions td at the target and sd at the
% source (a derivative with respect to y is minus one with respect to x)
Dg = @(td,sd) (-1)^numel(sd)*chnk.radpp.contract(comps,dim,[td sd]);

submat = zeros(info.opdims(1)*nt,info.opdims(2)*ns);

switch info.type

    case 's'
        submat = Dg({},{});

    case 'd'
        sn = sdir(normals(srcinfo),dim);
        submat = Dg({},{sn});

    case 'sp'
        tn = tdir(normals(targinfo),dim);
        submat = Dg({tn},{});

    case 'dp'
        tn = tdir(normals(targinfo),dim);
        sn = sdir(normals(srcinfo),dim);
        submat = Dg({tn},{sn});

    case 'c'
        sn = sdir(normals(srcinfo),dim);
        submat = coef(1)*Dg({},{sn}) + coef(2)*Dg({},{});

    case 'cp'
        tn = tdir(normals(targinfo),dim);
        sn = sdir(normals(srcinfo),dim);
        submat = coef(1)*Dg({tn},{sn}) + coef(2)*Dg({tn},{});

    case 'sg'
        for i = 1:dim
            submat(i:dim:end,:) = Dg({unitdir(i,dim)},{});
        end

    case 'dg'
        sn = sdir(normals(srcinfo),dim);
        for i = 1:dim
            submat(i:dim:end,:) = Dg({unitdir(i,dim)},{sn});
        end

    case 'cg'
        sn = sdir(normals(srcinfo),dim);
        for i = 1:dim
            submat(i:dim:end,:) = ...
                coef(1)*Dg({unitdir(i,dim)},{sn}) + ...
                coef(2)*Dg({unitdir(i,dim)},{});
        end

    case 'st'
        % tangent from the normal, as in chnk.lap2d.kern
        tn = normals(targinfo);
        tt = tdir([-tn(2,:); tn(1,:)],dim);
        submat = Dg({tt},{});

    case 'hilb'
        % two times the adjoint of st, as in chnk.lap2d.kern
        sn = sdir(normals(srcinfo),dim);
        submat = 2*(Dg({unitdir(1,dim)},{}).*sn{2} ...
            - Dg({unitdir(2,dim)},{}).*sn{1});

    case {'s2trans', 'clamped_bc'}
        tn = tdir(normals(targinfo),dim);
        submat(1:2:end,:) = Dg({},{});
        submat(2:2:end,:) = Dg({tn},{});

    case 'd2trans'
        tn = tdir(normals(targinfo),dim);
        sn = sdir(normals(srcinfo),dim);
        submat(1:2:end,:) = Dg({},{sn});
        submat(2:2:end,:) = Dg({tn},{sn});

    case 'c2trans'
        tn = tdir(normals(targinfo),dim);
        sn = sdir(normals(srcinfo),dim);
        submat(1:2:end,:) = coef(1,1)*Dg({},{sn}) + coef(1,2)*Dg({},{});
        submat(2:2:end,:) = coef(2,1)*Dg({tn},{sn}) + coef(2,2)*Dg({tn},{});

    case 'trans'
        tn = tdir(normals(targinfo),dim);
        sn = sdir(normals(srcinfo),dim);
        submat(1:2:end,1:2:end) = coef(1,1)*Dg({},{sn});
        submat(1:2:end,2:2:end) = coef(1,2)*Dg({},{});
        submat(2:2:end,1:2:end) = coef(2,1)*Dg({tn},{sn});
        submat(2:2:end,2:2:end) = coef(2,2)*Dg({tn},{});

    case 'trans_rep'
        sn = sdir(normals(srcinfo),dim);
        submat(:,1:2:end) = coef(1)*Dg({},{sn});
        submat(:,2:2:end) = coef(2)*Dg({},{});

    case 'trans_rep_prime'
        tn = tdir(normals(targinfo),dim);
        sn = sdir(normals(srcinfo),dim);
        submat(:,1:2:end) = coef(1)*Dg({tn},{sn});
        submat(:,2:2:end) = coef(2)*Dg({tn},{});

    case 'trans_rep_grad'
        sn = sdir(normals(srcinfo),dim);
        for i = 1:dim
            submat(i:dim:end,1:2:end) = coef(1)*Dg({unitdir(i,dim)},{sn});
            submat(i:dim:end,2:2:end) = coef(2)*Dg({unitdir(i,dim)},{});
        end

    case 'clamped_eval'
        sn = sdir(normals(srcinfo),dim);
        st = sdir(tangents(srcinfo),dim);
        submat(:,1:2:end) = Dg({},{sn,sn,sn}) + 3*Dg({},{sn,st,st});
        submat(:,2:2:end) = -Dg({},{sn,sn}) + Dg({},{st,st});

    case 'supported_bc'
        tn = tdir(normals(targinfo),dim);
        tt = tdir(tangents(targinfo),dim);
        submat(1:2:end,:) = Dg({},{});
        submat(2:2:end,:) = Dg({tn,tn},{}) + nu*Dg({tt,tt},{});

    case 'supported_eval'
        sn = sdir(normals(srcinfo),dim);
        st = sdir(tangents(srcinfo),dim);
        kappa = curvature(srcinfo);
        kp = srcinfo.data(1,:);
        a1 = 2 - nu;
        a2 = (nu-1)*(7+nu)/(3-nu);
        a3 = (1-nu)*(3+nu)/(1+nu);
        submat(:,1:2:end) = Dg({},{sn,sn,sn}) + a1*Dg({},{sn,st,st}) + ...
            a2*kappa.*Dg({},{sn,sn}) + a3*kp.*Dg({},{st});
        submat(:,2:2:end) = Dg({},{sn});

    case 'free_bc'
        tn = tdir(normals(targinfo),dim);
        tt = tdir(tangents(targinfo),dim);
        kappa = curvature(targinfo).';
        hnn = Dg({tn,tn},{});
        htt = Dg({tt,tt},{});
        submat(1:2:end,:) = hnn + nu*htt;
        submat(2:2:end,:) = Dg({tn,tn,tn},{}) + (2-nu)*Dg({tn,tt,tt},{}) ...
            + (1-nu)*kappa.*(htt - hnn);

    case 'custom'
        % user weights on the derivatives, entries in column-major order
        td2 = @(c) cellfun(@(a) tdir(a,dim),c,'UniformOutput',false);
        sd2 = @(c) cellfun(@(a) sdir(a,dim),c,'UniformOutput',false);
        D = @(td,sd) Dg(td2(td),sd2(sd));
        W = type.fun(D,srcinfo,targinfo);
        m = info.opdims(1);  n = info.opdims(2);
        for e = 1:m*n
            submat(mod(e-1,m)+1:m:end,ceil(e/m):n:end) = W(:,:,e);
        end

    case 'free_eval'
        sn = sdir(normals(srcinfo),dim);
        st = sdir(tangents(srcinfo),dim);
        submat(:,1:3:end) = Dg({},{sn});
        submat(:,2:3:end) = (1+nu)/2*Dg({},{st});
        submat(:,3:3:end) = Dg({},{});

end

end


function n = normals(p)
% unit normals of the points, from p.n or else p.d
if isfield(p,'n')
    n = p.n(:,:);
else
    n = chnk.normal2d(p);
end
end


function t = tangents(p)
% unit tangents of the points, from p.d
d = p.d(:,:);
t = d./sqrt(sum(d.^2,1));
end


function kappa = curvature(p)
% signed curvature (1,np) of the points, from p.d and p.d2
d = p.d(:,:);
d2 = p.d2(:,:);
kappa = (d(1,:).*d2(2,:) - d2(1,:).*d(2,:))./sqrt(sum(d.^2,1)).^3;
end


function v = tdir(a,dim)
% direction components of a vector field at the targets, a cell of dim
% (nt,1) arrays
v = cell(1,dim);
for i = 1:dim
    if i <= size(a,1), v{i} = a(i,:).'; else, v{i} = 0; end
end
end


function v = sdir(a,dim)
% direction components of a vector field at the sources, a cell of dim
% (1,ns) arrays
v = cell(1,dim);
for i = 1:dim
    if i <= size(a,1), v{i} = a(i,:); else, v{i} = 0; end
end
end


function v = unitdir(i,dim)
% the i-th coordinate direction
v = num2cell(double((1:dim) == i));
end
