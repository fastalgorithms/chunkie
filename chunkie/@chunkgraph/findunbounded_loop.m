function [loops_in,loops_out] = findunbounded_loop(cgrph,loops)
%FINDUNBOUNDED after having computed all the loops, loops, of a 
% chunkgraph, cgrph, findunbounded identifies `exterior' loops (loop_out) 
% and `interior' loops (loop_in).
%
% Syntax: [rgn] = findunbounded_loop(cgrph,rgn);
%
% Input:
%   cgrph  - chunkgraph or chunkgraph_per object
%   loops  - the loop cell array containing edge indices of each loop
%
% Output:
%   loop_in, loop_out.
%  
%
%

% author: Jeremy Hoskins

    iunbound = 1;
    
    inds_out = [];

    for ii=1:numel(loops)
        
        
        
        edges = loops{ii};
        theta = 0;
        rends = zeros([4,numel(edges)]);
        tends = zeros([4,numel(edges)]);
        
        for jj =1:numel(edges)
            echnk = cgrph.echnks(abs(edges(jj)));     
            [allrends,alltends] = chunkends(echnk);
            ts1 = alltends(:,1,:);
            ts2 = alltends(:,2,:);
            angs1 = atan2(real(ts1(2,:)),real(ts1(1,:)));
            angs2 = atan2(real(ts2(2,:)),real(ts2(1,:)));
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
            ang1 = atan2(real(tv1(2)),real(tv1(1)));
            ang2 = atan2(real(tv2(2)),real(tv2(1)));
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
            inds_out = [inds_out,ii];
        end
        
    end
    loops_out = loops(inds_out);
    loops_in  = loops;
    loops_in(inds_out) = [];
end

