function [np,alpha1,alpha2] = get_alphas_np_gui(chnkr,clmparams)

    if isfield(clmparams,'k')
      k = clmparams.k;
    end
    if isfield(clmparams,'c')
      c = clmparams.c;
    end
    if isfield(clmparams,'coef')
      coef = clmparams.coef;
    end
    if isfield(clmparams,'cpars')
      cpars = clmparams.cpars;
    end

    if isfield(clmparams,'ncorner')
      ncorner = clmparams.ncorner;
    end
    if isfield(clmparams,'corners')
      corners = clmparams.corners;
    end
    if isfield(clmparams,'issymmetric')
      issymmetric = clmparams.issymmetric;
    end

    if isfield(clmparams, 'nch')
      nch = clmparams.nch;
    end

    ncurve = length(chnkr);


% total number of discretization points
    ngl = chnkr(1).k;
    np = sum(nch(1:ncurve))*ngl;


% diagonal constant for each curve
    alpha1 = zeros(1,ncurve);
    alpha2 = zeros(1,ncurve);
    for i=1:ncurve
      alpha1(i) = 2/(coef(c(1,i))+coef(c(2,i)));
      alpha2(i) = 2/(1/coef(c(1,i))+1/coef(c(2,i)));
    end
    
end
