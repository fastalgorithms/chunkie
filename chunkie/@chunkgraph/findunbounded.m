function [rgn] = findunbounded(cgrph,rgn)
%FINDUNBOUNDED after having computed all the regions, rgn, of a 
% chunkgraph, cgrph, findunbounded identifies an unbounded region and 
% moves it to the start of the rgn cell array.
% NOTE: this routine is designed for finding the complement of 
% simply-connected regions. For non-simply connected, two rgn arrays 
% should be built for each component and findunbounded should be called 
% on both separately. The two can then be merged manually.
%
% Syntax: [rgn] = findunbounded(cgrph,rgn);
%
% Input:
%   cgrph  - chunkgraph object
%   rgn    - the rgn cell array containing edge indices of each region
%
% Output:
%   rgn    - the same cell array as rgn but with the unbounded region 
%            in the first entry.
%  
%
%

% author: Jeremy Hoskins

    iunbound = 1;
    
    for ii=1:numel(rgn)
        
        
        
        edges = rgn{ii}{1};
        theta = 0;
        rends = zeros([4,numel(edges)]);
        tends = zeros([4,numel(edges)]);
        
        for jj =1:numel(edges)
            echnk = cgrph.echnks(abs(edges(jj)));     
            [allrends,alltends] = chunkends(echnk);
            ts1 = alltends(:,1,:);
            ts2 = alltends(:,2,:);
            angs1 = atan2(ts1(2,:),ts1(1,:));
            angs2 = atan2(ts2(2,:),ts2(1,:));
            angdiffs = angs2-angs1;
            angdiffs(angdiffs>pi) = angdiffs(angdiffs>pi)-2*pi;
            angdiffs(angdiffs<-pi) = angdiffs(angdiffs<-pi)+2*pi;
            tchnk = sign(edges(jj))*sum(angdiffs);
            theta = theta + tchnk;
            
           	[rend,tend] = chunkends(echnk,1);
          	rend1 = rend(:,1);
          	tend1 = tend(:,1);
                
           	[rend,tend] = chunkends(echnk,echnk.nch);
            rend2 = rend(:,2);
           	tend2 = tend(:,2);
            
            if (edges(jj)>0)
            rends(:,jj) = [rend1;rend2];
            tends(:,jj) = [tend1;tend2];
            else
          	rends(:,jj) = [rend2;rend1];
            tends(:,jj) = [-tend2;-tend1];    
            end
        end
        
        tends = [tends,tends(:,1)];
        
        angsum = 0;
        for jj=1:numel(edges)
            tv1 = tends(3:4,jj);
            tv2 = tends(1:2,jj+1);
            ang1 = atan2(tv1(2),tv1(1));
            ang2 = atan2(tv2(2),tv2(1));
            angdiff = ang2-ang1;
            if (angdiff < -pi)
                angdiff = angdiff + 2*pi;
            end    
            if (angdiff >= pi)
                angdiff = angdiff - 2*pi;
            end
            angsum = angsum + angdiff;
            
        end  
        
        tot_ang = angsum + theta;
        if (tot_ang>pi)
            iunbound = ii;
        end
        
    end
    
    rgn([1,iunbound]) = rgn([iunbound,1]);
end

