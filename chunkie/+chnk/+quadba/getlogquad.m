function [xs,ws,nskip] = getlogquad(norder)

filename= '+chnk/+quadba/logwhts.mat';

persistent whts_dict
if isempty(whts_dict)
    load(filename,'whts_dict');
    assert(~isempty(whts_dict),'failed to load!')
end

assert(isKey(whts_dict,norder),'not an available integration order');

whts = whts_dict(norder);

xs = whts.xs; ws = whts.ws; nskip = whts.nskip;

end
