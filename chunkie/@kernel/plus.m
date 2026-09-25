function f = plus(f,g)
% + Pointwise addition for kernel class
%
% Currently only supported for adding two kernel class objects.

if (isa(g,'kernel') && isa(f,'kernel'))
  assert(f.opdims(1) == g.opdims(1) && f.opdims(2) == g.opdims(2), ...
      'kernel dimensions must agree to add');

  % One part per singularity type, so chunkermat can use a matching rule for each
  if ~isempty(f.sing) && ~isempty(g.sing)
    fparts = f.parts;
    gparts = g.parts;
    if isempty(fparts), fparts = struct(lower(f.sing), f); end
    if isempty(gparts), gparts = struct(lower(g.sing), g); end
    if f.iszero, fparts = struct(); end
    if g.iszero, gparts = struct(); end

    sings = union(fieldnames(fparts), fieldnames(gparts));
    if any(strcmp(sings, 'smooth')),    f.sing = 'smooth';    end
    if any(strcmp(sings, 'removable')), f.sing = 'removable'; end
    if any(strcmp(sings, 'log')),       f.sing = 'log';       end
    if any(strcmp(sings, 'pv')),        f.sing = 'pv';        end
    if any(strcmp(sings, 'hs')),        f.sing = 'hs';        end

    f.parts = [];
    if numel(sings) > 1
      for sname = fieldnames(gparts)'
        s = sname{1};
        if isfield(fparts, s)
          fparts.(s) = fparts.(s) + gparts.(s);
        else
          fparts.(s) = gparts.(s);
        end
      end
      f.parts = fparts;
    end
  else
    f.parts = [];
  end

  f.name = ['custom ',f.name,' ',g.name];

  if(isa(f.shifted_eval, 'function_handle'))
    if(isa(g.shifted_eval, 'function_handle'))
        f.shifted_eval = @(varargin) f.shifted_eval(varargin{:})+g.shifted_eval(varargin{:});
    else
        f.shifted_eval = @(varargin) f.shifted_eval(varargin{:})+g.eval(varargin{1:2});
    end
  else
    if(isa(g.shifted_eval, 'function_handle'))
        f.shifted_eval = @(varargin) f.eval(varargin{1:2})+g.shifted_eval(varargin{:});
    else
        f.shifted_eval = [];
    end
  end

  f.splitinfo = kernel.combine_splitinfo(f.splitinfo, g.splitinfo, 1, 1);

  f.eval = @(varargin) g.eval(varargin{:}) + f.eval(varargin{:});
  if (isa(g.fmm,'function_handle') && isa(f.fmm,'function_handle'))
    f.fmm = @(varargin) g.fmm(varargin{:}) + f.fmm(varargin{:});
  else
    f.fmm = [];
  end

  if or(f.isnan,g.isnan)
      f = kernel.nans(f.opdims(1),f.opdims(2));
  end
  if and(f.iszero,g.iszero)
      f.iszero = true;
  else
      f.iszero = false;
  end
else
    error('KERNEL:plus:invalid', ...
       'F and G must be kernel class objects');
end
end

