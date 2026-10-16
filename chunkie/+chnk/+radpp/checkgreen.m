function checkgreen(green,nderiv,dim,name)
%CHNK.RADPP.CHECKGREEN check that a Green's function is radial and has
% consistent derivatives
%
% Syntax: chnk.radpp.checkgreen(green,nderiv,dim,name)
%
% Errors unless green returns nderiv + 1 outputs of the right sizes, the
% outputs are the derivative tensors of a radial function, and each
% output matches finite differences of the one before. name is the kernel
% type, for the error messages.
%
% see also KERNEL.RADGREEN, CHNK.RADPP.PROFILES

err0 = 'KERNEL.RADGREEN: ';
names = {'val', 'grad', 'hess', 'third', 'fourth'};
need = strjoin(names(1:nderiv+1),', ');

nt = 3;  ns = 2;
s = randn(dim,ns);  t = randn(dim,nt) + [2; zeros(dim-1,1)];
out = cell(1,nderiv+1);
try
    [out{:}] = green(s,t);
catch ME
    error([err0 'type ''%s'' needs [%s] from green: %s'],name,need, ...
        ME.message);
end

if ~isequal(size(out{1}),[nt ns])
    error([err0 'val of green should be %s, not %s.'], ...
        mat2str([nt ns]),mat2str(size(out{1})));
end
C = out{1}(:);
for k = 1:nderiv
    try
        c = chnk.radpp.unpack(out{k+1},dim,k);
    catch ME
        error([err0 '%s of green: %s'],names{k+1},ME.message);
    end
    if size(c,1) ~= nt*ns
        error([err0 '%s of green should have %d points, not %d.'], ...
            names{k+1},nt*ns,size(c,1));
    end
    C = [C, c];
end

% radial: the tensors are those of the on-axis profiles
dx = t(1,:).' - s(1,:);
dy = t(2,:).' - s(2,:);
dz = 0*dx;
if dim == 3, dz = t(3,:).' - s(3,:); end
r = sqrt(dx.^2 + dy.^2 + dz.^2);
prof = chnk.radpp.profiles(green,r(:),nderiv,dim);
P = chnk.radpp.rotate(dx,dy,dz,prof,nderiv,dim);
if norm(P - C,'fro') > 1e-8*norm(C,'fro')
    error([err0 'green is not radial, or its derivatives are not those ' ...
        'of an isotropic tensor.']);
end

% each derivative against finite differences of the one below
if nderiv == 0, return; end
y0 = [0.2; -0.1; 0.4];  y = y0(1:dim);
x0 = [0.5; 0.6; -0.3];  x = y + x0(1:dim);
h = 1e-5;
base = allcomps(green,y,x,nderiv,dim);
Cp = cell(1,dim);  Cm = cell(1,dim);
for j = 1:dim
    e = zeros(dim,1);  e(j) = h;
    Cp{j} = allcomps(green,y,x+e,nderiv-1,dim);
    Cm{j} = allcomps(green,y,x-e,nderiv-1,dim);
end

for k = 0:nderiv-1
    tup = chnk.radpp.tuples(dim,k);
    tup1 = chnk.radpp.tuples(dim,k+1);
    nk = size(tup,1);
    A = zeros(dim,nk);
    F = zeros(dim,nk);
    for j = 1:dim
        F(j,:) = (Cp{j}{k+1} - Cm{j}{k+1})/(2*h);
        for m = 1:nk
            pos = find(all(tup1 == sort([tup(m,:) j]),2));
            A(j,m) = base{k+2}(pos);
        end
    end
    rel = norm(A - F,'fro')/max(norm(A,'fro'),realmin);
    if rel > 1e-6
        error([err0 '%s of green does not match finite differences ' ...
            'of %s (%.1e).'],names{k+2},names{k+1},rel);
    end
end

end


function C = allcomps(green,y,x,nderiv,dim)
% value and derivative components of green at one source and target,
% C{k+1} is (1,ncomp) for the derivatives of order k
out = cell(1,nderiv+1);
[out{:}] = green(y,x);
C = cell(1,nderiv+1);
C{1} = out{1}(:).';
for k = 1:nderiv
    C{k+1} = chnk.radpp.unpack(out{k+1},dim,k);
end
end
