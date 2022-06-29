function plot(obj,varargin)


ifhold = ishold();

echnks =  obj.echnks;
regions = obj.regions;

hold on

for ii=1:numel(regions)
   
    for jj=2:numel(regions{ii})

        
        rs = [];
        for ijk=1:numel(regions{ii}{jj}{1})
            enum = regions{ii}{jj}{1}(ijk);
            rchnk = echnks(abs(enum)).r;
            rchnk =rchnk(1:2,:);
            if (enum<0)
               rchnk = fliplr(rchnk); 
            end    
            rs = [rs,rchnk];
        end  
        plyrgn = polyshape(rs');
        
        for kk=2:numel(regions{ii}{jj})
        rs = [];
        for ijk=1:numel(regions{ii}{jj}{kk})
            enum = regions{ii}{jj}{kk}(ijk);
            rchnk = echnks(abs(enum)).r;
            rchnk =rchnk(1:2,:);
            if (enum<0)
               rchnk = fliplr(rchnk); 
            end    
            rs = [rs,rchnk];
        end
        
        plyrgnsub = polyshape(rs');
        plyrgn = subtract(plyrgn,plyrgnsub);
        end  
        plot(plyrgn);
    end    
    
end    

hold off

if ifhold
    hold on
end