% GENPATH_EX  Like GENPATH but excluding hidden and other special directories.
%
%    See also GENPATH.

function p = genpath_ex(root)
  p = '';
  files = dir(root);
  for i = 1:length(files)
    if ~files(i).isdir, continue; end
    d = files(i).name;
    if d(1) == '.' || strcmpi(d,'devtools') || strcmpi(d,'docs') || d(1) == '+' || d(1) == '@' || strcmpi(d,'clmtest') || strcmpi(d,'demo') , continue; end
    p = [p fullfile(root,d) pathsep];
    p = [p genpath_ex(fullfile(root,d))];
  end
end
