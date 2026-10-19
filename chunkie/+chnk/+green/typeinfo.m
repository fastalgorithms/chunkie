function info = typeinfo(type,dim)
%CHNK.GREEN.TYPEINFO name, operator dimensions and derivative order of the
% kernel types of KERNEL.GREEN and KERNEL.RADPP
%
% Syntax: info = chnk.green.typeinfo(type,dim)
%
% info.type   - the name of the type, aliases are resolved
% info.opdims - operator dimensions [m n] of the kernel
% info.nderiv - highest derivative of the Green's function needed, source
%               derivatives included
% info.param  - the parameter taken after the type: '' (none), 'coefs'
%               (2 coefficients), 'coefs22' (2x2 coefficients) or 'nu'
%
% see also KERNEL.GREEN, CHNK.GREEN.KERN

if nargin < 2 || isempty(dim), dim = 2; end

planar = false;

switch lower(type)
    case {'s', 'single'}
        [name,opdims,nderiv,param] = deal('s',[1 1],0,'');
    case {'d', 'double'}
        [name,opdims,nderiv,param] = deal('d',[1 1],1,'');
    case {'sp', 'sprime'}
        [name,opdims,nderiv,param] = deal('sp',[1 1],1,'');
    case {'dp', 'dprime'}
        [name,opdims,nderiv,param] = deal('dp',[1 1],2,'');
    case {'c', 'combined'}
        [name,opdims,nderiv,param] = deal('c',[1 1],1,'coefs');
    case {'cp', 'cprime'}
        [name,opdims,nderiv,param] = deal('cp',[1 1],2,'coefs');
    case {'sg', 'sgrad'}
        [name,opdims,nderiv,param] = deal('sg',[dim 1],1,'');
    case {'dg', 'dgrad'}
        [name,opdims,nderiv,param] = deal('dg',[dim 1],2,'');
    case {'cg', 'cgrad'}
        [name,opdims,nderiv,param] = deal('cg',[dim 1],2,'coefs');
    case {'st', 'stau'}
        [name,opdims,nderiv,param] = deal('st',[1 1],1,'');
        planar = true;
    case {'hilb', 'hilbert'}
        [name,opdims,nderiv,param] = deal('hilb',[1 1],1,'');
        planar = true;
    case {'s2trans', 's2t'}
        [name,opdims,nderiv,param] = deal('s2trans',[2 1],1,'');
    case {'d2trans', 'd2t'}
        [name,opdims,nderiv,param] = deal('d2trans',[2 1],2,'');
    case {'c2trans', 'c2t', 'c2tr'}
        [name,opdims,nderiv,param] = deal('c2trans',[2 1],2,'coefs22');
    case {'trans', 'all', 'trans_sys', 'tsys'}
        [name,opdims,nderiv,param] = deal('trans',[2 2],2,'coefs22');
    case {'trans_rep', 'trep'}
        [name,opdims,nderiv,param] = deal('trans_rep',[1 2],1,'coefs');
    case {'trans_rep_prime', 'trep_p', 'trans_rep_p'}
        [name,opdims,nderiv,param] = deal('trans_rep_prime',[1 2],2,'coefs');
    case {'trans_rep_grad', 'trep_g', 'trans_rep_g'}
        [name,opdims,nderiv,param] = deal('trans_rep_grad',[dim 2],2,'coefs');
    case {'clamped_bc', 'clamped', 'clamped_plate_bcs', 'clamped_plate_bc'}
        [name,opdims,nderiv,param] = deal('clamped_bc',[2 1],1,'');
    case {'clamped_eval', 'clamped_plate_eval'}
        [name,opdims,nderiv,param] = deal('clamped_eval',[1 2],3,'');
        planar = true;
    case {'supported_bc', 'supported', 'supported_plate_bcs', 'supported_plate_bc'}
        [name,opdims,nderiv,param] = deal('supported_bc',[2 1],2,'nu');
        planar = true;
    case {'supported_eval', 'supported_plate_eval'}
        [name,opdims,nderiv,param] = deal('supported_eval',[1 2],3,'nu');
        planar = true;
    case {'free_bc', 'free', 'free_plate_bcs', 'free_plate_bc'}
        [name,opdims,nderiv,param] = deal('free_bc',[2 1],3,'nu');
        planar = true;
    case {'free_eval', 'free_plate_eval'}
        [name,opdims,nderiv,param] = deal('free_eval',[1 3],1,'nu');
        planar = true;
    otherwise
        error('CHNK.GREEN.TYPEINFO: unknown kernel type ''%s''.',type);
end

if planar && dim ~= 2
    error('CHNK.GREEN.TYPEINFO: type ''%s'' needs dim = 2.',name);
end

info = struct('type',name,'opdims',opdims,'nderiv',nderiv,'param',param);

end
