function [obj] = balance(obj)
%BALANCE Given a vertex, with incident edges E1,E2,...,Ek, balance 
% makes sure that the arclength of the chunks of E1,...,EK, nearest the 
% vertex are within a factor of two of eachother.
%
% Syntax: [cgrph] = balance(cgrph);
%
% Input:
%   cgrph  - chunkgraph object
%
% Output:
%   cgrph  - chunkgraph object with refined chunks near vertices.
%  
%
%

% author: Jeremy Hoskins

    if (isfield(obj,'vstruc'))
        vstruc = obj.vstruc;
    else
        vstruc = procverts(obj);
    end
    
    echnks = obj.echnks;
    
    nverts = size(obj.verts,2);
    % Loop over all vertices in the cgraph structure
    for iii=1:nverts
        if isempty(vstruc{iii}{1})
            % if vertex is not connected, do nothing
            continue
        end
    vedge = vstruc{iii}{1};
    vsign = vstruc{iii}{2};
    
    ifdone = false;
    while (~ifdone)
    % find the arclengths of each panel on the edge incident to 
    % the vertex 
    
    parcl = zeros([numel(vedge),1]);
    pinds = zeros([numel(vedge),1]);
    
    for ii=1:numel(vedge)
        if (sign(vsign(ii)) == -1)
            ds =  obj.echnks(vedge(ii)).d(:,:,1);
            pinds(ii) = 1;
        else
            ds =  obj.echnks(vedge(ii)).d(:,:,end);
            pinds(ii) = size(obj.echnks(vedge(ii)).r,3);
        end    
        
        k  =  obj.echnks(vedge(ii)).k;

        wleg = echnks(vedge(ii)).wstor;

        arc = sum(sqrt(ds(1,:).^2+ds(2,:).^2).*wleg');
        parcl(ii) = arc;
          
    end    
    % find the ratio of the max to min panels of all edges
    % incident at the vertex and refine the panel with the 
    % max panel length
    amin = min(parcl);
    [amax,ind] = max(parcl);
    ind = ind(1);
    n2ref = floor(log(amax/amin)/log(2));
    
    
    if (n2ref >0)
        chnkr = obj.echnks(vedge(ind));
        opts = [];
        opts.lvlrfac = 'a';
        for ilev = 1:n2ref
            opts.splitchunks=([pinds(ind)]);
            % If it is the last chunk, then update pinds to 
            % ensure that the last chunk is refined at 
            % subsequent levels
            if (pinds(ind)>1)
                pinds(ind) = pinds(ind)+1;
            end
            % resort chnkr in case both end points need
            % to be refined in the wrong order
            chnkr = refine(chnkr,opts);
            chnkr = sort(chnkr);
        end
        obj.echnks(vedge(ind)) = chnkr;
        
    else
        ifdone = true;
    end
    
    end
    end
    
    
end
