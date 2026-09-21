function s = scale_splitinfo(s0, c)
%KERNEL.SCALE_SPLITINFO   Scale kernel-split metadata by scalar C, for C*K.

if ~(isstruct(s0) && isfield(s0,'functions') && ~isempty(s0.functions) ...
        && isfield(s0,'type') && isfield(s0,'action'))
    s = [];
    return
end

s = s0;
f0 = s0.functions;
s.functions = @(src,targ) cellfun(@(x) c*x, f0(src,targ), 'UniformOutput', false);

end
