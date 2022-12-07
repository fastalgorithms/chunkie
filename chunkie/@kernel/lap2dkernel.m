function obj = lap2dkernel(type,opts)
		     % return the laplace kernel of the requested type
  obj = emptykernel();
  if (strcmpi(type,'d') || strcmpi(type,'double') || strcmpi(type,'double layer'))
    obj.eval = @(s,t) chnk.lap2d.kern('d',s,t);
    obj.bdrysing = 'smooth';
  end
end

