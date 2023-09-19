function mat = buildmat(chnkr,fkern,opts,opdims,iglist)
% Build system matrix for rcip purposes
% 
% 
   ncurve = length(chnkr);
   nch = zeros(1,ncurve);
   for i=1:ncurve
       nch(i) = chnkr(i).nch;
   end
   ngl = chnkr(1).k;
   np = sum(nch(1:ncurve))*ngl;
   mat = zeros(np*opdims(2));
   
   chnkinfo = cell(1,ncurve);
   
   for i=1:ncurve
       chnkinfo{i} = [];
       chnkinfo{i}.r = chnkr(i).r(:,:);
       chnkinfo{i}.d = chnkr(i).d(:,:);
       chnkinfo{i}.d2 = chnkr(i).d2(:,:);
       chnkinfo{i}.n = chnkr(i).n(:,:);
       wts = chnkr(i).wts;
       wts2 = repmat( (wts(:)).', opdims(2), 1);
       wts2 = ( wts2(:) ).';
       chnkinfo{i}.wts = wts2;
   end
   for i=1:ncurve
       indi1 = (1:(opdims(1)*nch(i)*ngl))+sum(nch(1:i-1))*ngl*opdims(1);
       indi2 = (1:(opdims(2)*nch(i)*ngl))+sum(nch(1:i-1))*ngl*opdims(2);
       for j=1:ncurve
           if j==i
               jlist = [];
               if ~isempty(iglist)
                   jlist = iglist(:,j);
               end
               mat(indi1,indi2) = chunkermat(chnkr(i),fkern,opts,jlist);               
           else
               indj = (1:(opdims(2)*nch(j)*ngl))+ sum(nch(1:(j-1)))*ngl*opdims(2);
               mat(indi1,indj) = fkern(chnkinfo{j},chnkinfo{i}).*chnkinfo{j}.wts;
           end
       end
   end
end
