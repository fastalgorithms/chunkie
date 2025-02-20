function obj = makedatarows(obj,nrows)
%MAKEDATAROWS add new data rows to chunkgraph object
%
%
% Syntax: obj = makedatarow(obj, nrows);
%
% Input:
%   obj    - a chunkgraph object
%   nrows  - number of rows to be added
%
% Output:
%   obj    - output chunkgraph object with additional data rows

%
%
    if (nrows > 0)
        datatemp = obj.data;
        datadimold = obj.datadim;
        nch = sum(horzcat(obj.echnks.nch));
        obj.data = zeros(datadimold+nrows,obj.k,nch);
        obj.data(1:datadimold,:) = datatemp(:,:);
        obj.hasdata = true;
    else
        if (nrows < 0)
            warning('attempted to add negative rows, doing nothing');
        end
    end
end
