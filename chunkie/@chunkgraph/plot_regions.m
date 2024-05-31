function plot_regions(obj,iflegend)
%PLOT_REGIONS plots regions of a chunkgraph in 2 dimensions
% All regions in the chunkgraph are plotted in a different color.
%
% Syntax: plot_regions(cgrph)
%
% Input: 
%   cgrph - chunkgraph object
%
% Optional input:
%   iflegend - boolean (default true), include a Legend showing the region
%   numbers
%
% Output:
%   none 
%

% author: Jeremy Hoskins

if nargin < 2 || isempty(iflegend)
    iflegend = true;
end


ifhold = ishold();

echnks =  obj.echnks;
regions = obj.regions;

hold on

nr = numel(regions);
legtext = cell(max(1,nr-1),1);

for ii=2:numel(regions)
    legtext{ii-1} = "region " + num2str(ii);
           rs = [];
        for ijk=1:numel(regions{ii}{1})
            enum = regions{ii}{1}(ijk);
            rchnk = echnks(abs(enum)).r;
            rchnk =rchnk(1:2,:);
            if (enum<0)
               rchnk = fliplr(rchnk); 
            end    
            rs = [rs,rchnk];
        end  
        plyrgn = polyshape(rs.','Simplify',false);

    for jj=2:numel(regions{ii})

        for kk=2:numel(regions{ii})
        rs = [];
        for ijk=1:numel(regions{ii}{kk})
            enum = regions{ii}{kk}(ijk);
            rchnk = echnks(abs(enum)).r;
            rchnk =rchnk(1:2,:);
            if (enum<0)
               rchnk = fliplr(rchnk); 
            end    
            rs = [rs,rchnk];
        end
        
        plyrgnsub = polyshape(rs','Simplify',false);
        plyrgn = subtract(plyrgn,plyrgnsub);
        end  
        
    end    
    plot(plyrgn);
end    

if nr > 1 && iflegend
    legend(legtext);
end


hold off

if ifhold
    hold on
end