function load_and_save_bremer_quad_info(filename)
%LOAD_AND_SAVE_BREMER_QUAD_INFO
%
% 
%
%  orders available = [16,20,24,30,40,60];

iords = [16,20,24,30,40,60];

whts_vals = cell(length(iords),1);

for i = 1:length(iords)
    whts = [];
    [xs1,whts1,xs0,whts0] = getquads(iords(i));
    whts.xs1 = xs1;
    whts.whts1 = whts1;
    whts.xs0 = xs0;
    whts.whts0 = whts0;
    whts_vals{i} = whts;
end

whts_dict = containers.Map(iords,whts_vals,'UniformValues',false);

save(filename,'whts_dict');