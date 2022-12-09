function [r,d,d2] = linefunc(t,v1,v2)
    sizet = size(t);
    ts = t(:).';
    r  = v1 + (v2-v1)*ts;
    d  = (v2-v1)*ones(size(ts));
    d2 = zeros(size(d));
    r = reshape(r,[2,sizet]);
    d = reshape(d,[2,sizet]);
    d2= reshape(d2,[2,sizet]);
end
