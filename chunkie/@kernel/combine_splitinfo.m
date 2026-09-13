function s = combine_splitinfo(sf, sg, cf, cg)
%KERNEL.COMBINE_SPLITINFO   Merge kernel-split metadata for CF*F + CG*G.

if ~is_split(sf) || ~is_split(sg)
    s = [];
    return
end

s = struct();
s.type   = [sf.type(:).',   sg.type(:).'];
s.action = [sf.action(:).', sg.action(:).'];

ff = sf.functions;
gf = sg.functions;
s.functions = @(src,targ) [ ...
    cellfun_row(@(x) cf*x, ff(src,targ)), ...
    cellfun_row(@(x) cg*x, gf(src,targ))];

end

function tf = is_split(s)
tf = isstruct(s) && isfield(s,'functions') && ~isempty(s.functions) ...
    && isfield(s,'type') && isfield(s,'action');
end

function c = cellfun_row(fun, c)
c = reshape(cellfun(fun, c, 'UniformOutput', false), 1, []);
end
