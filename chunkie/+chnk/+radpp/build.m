function obj = build(green,rmax,type,args,isfit)
%CHNK.RADPP.BUILD kernel object from a radial Green's function, for
% KERNEL.RADGREEN and KERNEL.RADPP
%
% Syntax: obj = chnk.radpp.build(green,rmax,type,args,isfit)
%
% args is a cell array of the arguments after the type, see KERNEL.RADGREEN.
% If isfit the Green's function is replaced by a piecewise polynomial fit
% on rmax, see CHNK.RADPP.FIT.
%
% see also KERNEL.RADGREEN, KERNEL.RADPP

if isnumeric(green) || ~isa(green,'function_handle')
    error('CHNK.RADPP.BUILD: Green''s function must be a function handle.');
end
if isempty(type), type = 's'; end

% a positional parameter after the type, as in kernel.helm2d: the combined
% layer coefficients, or nu for the plate kernels
npos = find(cellfun(@ischar,args),1) - 1;
if isempty(npos), npos = numel(args); end
if npos > 1
    error('CHNK.RADPP.BUILD: One positional parameter (coefs or nu) after the type.');
end
par = [];
if npos == 1, par = args{1}; end
args(1:npos) = [];

o = struct('dim',2,'coefs',[],'nu',[],'sing',[],'opts',struct(), ...
    'nderiv',[],'opdims',[1 1]);
if mod(numel(args),2) ~= 0
    error('CHNK.RADPP.BUILD: Options must be name-value pairs.');
end
for i = 1:2:numel(args)
    name = lower(args{i});
    if ~isfield(o,name)
        error('CHNK.RADPP.BUILD: Unknown option ''%s''.',args{i});
    end
    o.(name) = args{i+1};
end
dim = o.dim;
if ~any(dim == [2 3])
    error('CHNK.RADPP.BUILD: dim must be 2 or 3.');
end

if isa(type,'function_handle')
    % custom kernel, the weights on the derivatives are given by type
    if isempty(o.nderiv) || ~any(o.nderiv == 0:4)
        error('CHNK.RADPP.BUILD: custom kernels need nderiv, 0 to 4.');
    end
    info = struct('type','custom','opdims',o.opdims,'nderiv',o.nderiv, ...
        'param','');
    ktype = struct('type','custom','fun',type,'nderiv',o.nderiv, ...
        'opdims',o.opdims);
else
    info = chnk.radpp.typeinfo(type,dim);
    ktype = info.type;
end

% the parameter of the type
switch info.param
    case ''
        if ~isempty(par) || ~isempty(o.coefs) || ~isempty(o.nu)
            error('CHNK.RADPP.BUILD: Type ''%s'' takes no coefs or nu.',info.type);
        end
        param = [];
    case 'nu'
        if ~isempty(par), o.nu = par; end
        if isempty(o.nu)
            error('CHNK.RADPP.BUILD: Type ''%s'' needs nu.',info.type);
        end
        param = o.nu;
    otherwise
        if ~isempty(par), o.coefs = par; end
        if isempty(o.coefs)
            if strcmp(info.param,'coefs')
                warning('Missing coefficients. Defaulting to [1 1].');
                o.coefs = [1 1];
            else
                warning('Missing coefficients. Defaulting to [1,1;1,1].');
                o.coefs = ones(2);
            end
        end
        param = o.coefs;
end

nfit = info.nderiv;
exactprof = @(r) chnk.radpp.profiles(green,r,nfit,dim);

obj = kernel();
obj.opdims = info.opdims;
obj.type = info.type;
obj.params.dim = dim;
obj.params.green = green;
obj.params.param = param;

if isfit
    % fit the profiles on [rmin rmax]
    rmin = 1e-12;
    if numel(rmax) == 2
        rmin = rmax(1);
        rmax = rmax(2);
    end
    fitopts = o.opts;
    if ~isfield(fitopts,'eps'), fitopts.eps = 1e-13; end
    fit = chnk.radpp.fit(exactprof,rmin,rmax,fitopts);
    prof = @(r) chnk.radpp.fiteval(r,fit.breaks,fit.coefs);

    obj.name = 'radpp';
    obj.params.rmin = rmin;
    obj.params.rmax = rmax;
    obj.params.fit = fit;
else
    chnk.radpp.checkgreen(green,nfit,dim,info.type);
    prof = exactprof;

    obj.name = 'radgreen';
end

obj.eval = @(s,t) chnk.radpp.kern(prof,dim,s,t,ktype,param);
obj.fmm = [];

% singularity type from the exact kernel
if isempty(o.sing)
    exacteval = @(s,t) chnk.radpp.kern(exactprof,dim,s,t,ktype,param);
    obj.sing = chnk.radpp.singtype(exacteval,info.opdims,dim);
else
    obj.sing = o.sing;
end

end
