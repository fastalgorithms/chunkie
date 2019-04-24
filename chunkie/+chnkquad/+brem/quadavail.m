function iords = quadavail()

filename= '+quad/+brem/whts.mat';

persistent whts_dict
if isempty(whts_dict)
    load(filename,'whts_dict');
    assert(~isempty(whts_dict),'failed to load!')
end

iords = keys(whts_dict);

end