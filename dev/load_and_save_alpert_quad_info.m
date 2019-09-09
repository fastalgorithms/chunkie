function load_and_save_alpert_quad_info(filename)
%LOAD_AND_SAVE_ALPERT_QUAD_INFO
%
% 
%
%  orders available = [0,2,3,4,8,16];

iords = [0,2,3,4,8,16];

whts_vals = cell(length(iords),1);

for i = 1:length(iords)
    whts = [];
    [xs,ws,nskip] = getquads_ba(iords(i));
    whts.xs = xs;
    whts.ws = ws;
    whts.nskip = nskip;
    whts_vals{i} = whts;
end

whts_dict = containers.Map(iords,whts_vals,'UniformValues',false);

save(filename,'whts_dict');