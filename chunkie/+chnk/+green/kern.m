function submat = kern(green,dim,srcinfo,targinfo,type,varargin)
%CHNK.GREEN.KERN kernels of a translation invariant Green's function
%
% Syntax: submat = chnk.green.kern(green,dim,srcinfo,targinfo,type,varargin)
%
% G(x,y) = g(x-y) is given by green, see KERNEL.GREEN. Source derivatives
% are minus target derivatives. S = G, D = G_{n_y}, S' = G_{n_x} and
% D' = G_{n_x n_y}.
%
% Input:
%   green - function handle, the Green's function
%   dim - 2 or 3, the dimension of the points
%   srcinfo, targinfo - sources and targets in ptinfo struct format
%   type - string, determines kernel type
%                type == 's', 'd', 'sp', 'dp', S, D, S', D'
%                type == 'c', 'cp', coef(1) D + coef(2) S and its normal
%                       derivative
%                type == 'sg', 'dg', 'cg', the gradients of S, D and C
%                type == 'st', 'hilb', tangential derivative of S and two
%                       times its adjoint (dim = 2)
%                type == 's2trans', 'd2trans', 'c2trans', 'trans',
%                       'trans_rep', 'trans_rep_prime', 'trans_rep_grad',
%                       as in CHNK.HELM2D.KERN
%                type == 'clamped_bc', 'clamped_eval', 'supported_bc',
%                       'supported_eval', 'free_bc', 'free_eval', plate
%                       kernels as in CHNK.FLEX2D.KERN (dim = 2)
%         or a custom struct, fun(D,srcinfo,targinfo) with D(td,sd) being
%                the derivative of G along the directions td and sd
%   varargin{1} - coef, 2x2 for c2trans and trans, or nu for the plate
%                kernels

% Output:
%   submat - the evaluation of the selected kernel for the provided
%            sources and targets.
%
% see also KERNEL.GREEN, KERNEL.RADPP, CHNK.GREEN.TYPEINFO

if isstruct(type)
    info = struct('type','custom','opdims',type.opdims, ...
        'nderiv',type.nderiv,'param','');
else
    info = chnk.green.typeinfo(type,dim);
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

ns = size(srcinfo.r(:,:),2);
nt = size(targinfo.r(:,:),2);

% derivative tensors of G with respect to the target, (nt,ns,ncomp)
nd = info.nderiv;
out = cell(1,nd+1);
[out{:}] = green(srcinfo.r(:,:),targinfo.r(:,:));
comps = cell(1,nd+1);
comps{1} = reshape(out{1},nt,ns);
for k = 1:nd
    c = chnk.green.unpack(out{k+1},dim,k);
    comps{k+1} = reshape(c,nt,ns,size(c,2));
end

% derivative of G along the directions td at the target and sd at the
% source (a derivative with respect to y is minus one with respect to x)
Dg = @(td,sd) (-1)^numel(sd)*contract(comps,dim,[td sd]);

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


function val = contract(comps,dim,vs)
% the derivative tensor of order k = numel(vs) contracted with the
% directions vs{j} = {vx,vy[,vz]}, sum D_{i1...ik} v1_{i1} ... vk_{ik}

k = numel(vs);
c = comps{k+1};
if k == 0
    val = c;
    return
end
tup = nchoosek(1:dim+k-1,k) - (0:k-1);
val = 0;
for m = 1:dim^k
    idx = cell(1,k);
    [idx{:}] = ind2sub([repmat(dim,1,k) 1],m);
    idx = [idx{:}];
    term = c(:,:,all(tup == sort(idx),2));
    for j = 1:k
        term = term.*vs{j}{idx(j)};
    end
    val = val + term;
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
