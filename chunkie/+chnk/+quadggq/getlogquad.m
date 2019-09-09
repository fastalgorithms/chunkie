function [xs1,whts1,xs0,whts0] = getlogquad(k)

filename= '+chnk/+quadggq/logwhts.mat';

persistent whts_dict
if isempty(whts_dict)
    load(filename,'whts_dict');
    assert(~isempty(whts_dict),'failed to load!')
end

assert(isKey(whts_dict,k),'not an available integration order');

whts = whts_dict(k);

xs1 = whts.xs1; whts1 = whts.whts1;
xs0 = whts.xs0; whts0 = whts.whts0;

end
