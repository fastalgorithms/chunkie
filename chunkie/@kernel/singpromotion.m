function sing = singpromotion(sing1,sing2)
%SINGPROMOTION promote singularity type when combining kernels using 
% hierarchy  
%
% Behavior: 
% - the hierarchy is based on available GGQ rules. smooth < log < pv < hs 
% - default singularity is empty for kernels. empty gets promoted to
% whatever the other singularity is if it's non-empty. two non-empty result
% in empty. warns when combining empty with non-empty.
% - if string other than hierarchy and the strings match, that is
% preserved. if two non-matching strings with any not in hierarchy, returns
% 'unknown'. 

% if they match don't think
if strcmpi(sing1,sing2) || (isempty(sing1) && isempty(sing2))
    sing = sing1;
    return
end

% try to use hierarchy 
hierarchy = {'smooth','log','pv','hs'};

index1 = find(strcmpi(sing1,hierarchy),1);
index2 = find(strcmpi(sing2,hierarchy),1);
if ~isempty(index1) && ~isempty(index2)
    sing = hierarchy{max(index1,index2)};
    return
end

% promote empty singularities to the non-empty type.
if isempty(sing1)
    sing = sing2;
    warning('ADDING KERNELS: empty singularity promoted to %s\n',sing);
    return
elseif isempty(sing2)
    sing = sing1;
    warning('ADDING KERNELS: empty singularity promoted to %s\n',sing);
    return
end

% if reached here one singularity type (or both) is out of scope and they
% do not match 

warning(['ADDING KERNELS: unclear how to combine %s and %s' ...
    ' singularity types, returning unknown\n'],sing1,sing2);
sing = 'unknown';

end