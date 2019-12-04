function iords = logavail()

filename= '+chnk/+quadggq/logwhts.mat';

persistent whts_dict
if isempty(whts_dict)
    load(filename,'whts_dict');
    assert(~isempty(whts_dict),'failed to load!')
end

iords = cell2mat(keys(whts_dict));

end
