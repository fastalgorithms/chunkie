function plot_regions(obj,iflabel)
%PLOT_REGIONS plots regions of a chunkgraph in 2 dimensions
% All regions in the chunkgraph are plotted in a different color.
%
% Syntax: plot_regions(cgrph)
%
% Input: 
%   cgrph - chunkgraph object
%
% Optional input:
%   iflabel - integer (default 2), include a Legend showing the region
%     numbers and label edges if 2, include only the region legend if 1, 
%     no labels if 0.
%
% Output:
%   none 
%

% author: Jeremy Hoskins

if nargin < 2 || isempty(iflabel)
    iflabel = 2;
end


ifhold = ishold();

echnks =  obj.echnks;
regions = obj.regions;


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
    hold on
end    

if nr > 1 && iflabel > 0
    legend(legtext,'AutoUpdate','off');
end

plot(obj,'k-');

if iflabel > 1
    rmin = min(obj.r(:,:),[],2);
    diam = max(obj.r(:,:)-rmin,[],2);
    rmin = rmin-0.1*diam;
    text(rmin(1),rmin(2),'region 1');

    nedge = length(obj.echnks);
    for j = 1:nedge
        wts = obj.echnks(j).wts;
        l1 = sum(wts(:));
        l2 = cumsum(wts(:));
        ind = find(l2 > l1/2,true,'first');
        r = obj.echnks(j).r(:,ind);
        text(r(1),r(2),"edge "+num2str(j),'HorizontalAlignment','center')
    end
end

hold off

if ifhold
    hold on
end