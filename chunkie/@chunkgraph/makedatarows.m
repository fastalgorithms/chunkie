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
        for ii = 1:length(obj.echnks)
            obj.echnks(ii) = makedatarows(obj.echnks(ii),nrows);
        end
    else
        if (nrows < 0)
            warning('attempted to add negative rows, doing nothing');
        end
    end
end
