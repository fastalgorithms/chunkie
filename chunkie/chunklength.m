function [len,dsdt] = chunklength(fcurve,a,b,xs,ws)
    
    nout = 3;
    out = cell(nout,1);
    ts = a+(b-a)*(xs+1)/2;
    [out{:}] = fcurve(ts);
    dsdt = sqrt(sum(out{2}.^2,1));
    len = dot(dsdt,ws)*(b-a)/2;
 end

