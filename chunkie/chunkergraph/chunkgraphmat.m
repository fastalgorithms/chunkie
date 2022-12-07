function [sysmat] = chunkgraphmat(chnklist,fkernarray,opts,ilist)


    indscol = [];
    indsrow = [];
    
    % determine operator dimensions using first two points

    i = 1;
    indsrow(1) = 1;
    indscol(1) = 1;
    
    for j=1:numel(chnklist)
        
        chnkri = chnklist(i);
        chnkrj = chnklist(j);
     	srcinfo = []; targinfo = [];
        srcinfo.r = chnkrj.r(:,1); srcinfo.d = chnkrj.d(:,1); 
        srcinfo.d2 = chnkrj.d2(:,1); srcinfo.n = chnkrj.n(:,1);
      	targinfo.r = chnkri.r(:,2); targinfo.d = chnkri.d(:,2); 
    	targinfo.d2 = chnkri.d2(:,2); targinfo.n = chnkri.n(:,2);

     	ftemp = fkernarray(srcinfo,targinfo,i,j);
      	opdims = size(ftemp);
        
        indsrow(j+1) = indsrow(j) + opdims(2)*size(chnkrj.r(:,:),2);   
    end

    j = 1;

    for i=1:numel(chnklist)
        
        chnkri = chnklist(i);
        chnkrj = chnklist(j);
     	srcinfo = []; targinfo = [];
        srcinfo.r = chnkrj.r(:,1); srcinfo.d = chnkrj.d(:,1); 
        srcinfo.d2 = chnkrj.d2(:,1); srcinfo.n = chnkrj.n(:,1);
      	targinfo.r = chnkri.r(:,2); targinfo.d = chnkri.d(:,2); 
    	targinfo.d2 = chnkri.d2(:,2); targinfo.n = chnkri.n(:,2);

     	ftemp = fkernarray(srcinfo,targinfo,i,j);
      	opdims = size(ftemp);
        
        indscol(i+1) = indscol(i) + opdims(1)*size(chnkri.r(:,:),2);   
    end
    
    
    indsrow
    indscol
    %%% now generate kernel
    
    sysmat = zeros(indsrow(end)-1,indscol(end)-1);
    
    for i=1:numel(chnklist)
        for j=1:numel(chnklist)
            if (i ~= j)
                kern = @(s,t) fkernarray(s,t,i,j);
                chnkrj = chnklist(j);
                chnkri = chnklist(i);
                targinfo = [];
                targinfo.r = chnkri.r(:,:); targinfo.d = chnkri.d(:,:); 
                targinfo.d2 = chnkri.d2(:,:); targinfo.n = chnkri.n(:,:);
                opdims
                
                mat = chunkerkernevalmat(chnkrj,kern,targinfo,opts);  
                [n1,n2] = size(mat);
                i1 = indsrow(i);
                i2 = indscol(j);
                sysmat(i1+(1:n1),i2+(1:n2)) = mat;
            end   
        end
    end
    
            
    
end

