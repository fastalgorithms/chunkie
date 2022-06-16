function [M,RG] = get_mat_gui_clm_cases(chnkr,clmparams,icase,opts)

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

    
    
    
    % number of Gauss-Legendre nodes on each chunk
    ngl = 16;
    [glnodes,glwts] = lege.exps(ngl);


    % Treat the representation as if it were a 2x2 operator so that four layer 
    % potentials D, S, D', S' can be evaluated together. This will avoid 
    % redundant evaluation of Hankel functions.
    opdims(1)=2;opdims(2)=2;

    % set up GGQ machinery
    [logquad] = chnk.quadggq.setuplogquad(ngl,opdims);

    % log correction for self interaction in kernel-split style
    isclosed = 0;
    hlocal = [1];
    LogC = chnk.quadjh.setuplogquad(hlocal,ngl,isclosed,opdims);
    logquad.LogC  = LogC;

    isrcip = 1;

    % build the system matrix
    rpars = [];
    rpars.k = k;
    rpars.c = c;
    rpars.coef = coef;

    disp(' ')
    disp('Step 1: build the system matrix directly.')
    start = tic;
    ilist=[];
    [M,np,~,~] = clm.buildmat_fast(chnkr,rpars,opts,opdims,glwts,ilist,logquad);
    if ~isrcip, M = M + eye(2*np); end
    dt = toc(start);

    disp(['System matrix construction time = ', num2str(dt), ' seconds'])
    %fprintf('%5.2f seconds : time to assemble matrix\n',t1)
    
    
    nonsmoothonly = false;  
    if isfield(opts,'nonsmoothonly')
      nonsmoothonly = opts.nonsmoothonly;
    end

    % compute the preconditioner R in the RCIP method for triple junctions

    if isrcip && ncorner>0
      disp(' ')
      disp('Step 2: compute the preconditioners for corners.')
      start = tic;

      if isreal(k)
        opts.quad = 'jhlog';
        hlocal1 = [0.5, 0.5, 1];
        hlocal0 = [1, 0.5, 0.5];

        LogC0 = chnk.quadjh.setuplogquad(hlocal0,ngl,isclosed,opdims);
        LogC1 = chnk.quadjh.setuplogquad(hlocal1,ngl,isclosed,opdims);

        logquad.LogC0 = LogC0;
        logquad.LogC1 = LogC1;
      end

      inds = [0, cumsum(nch)];

      ndim = opdims(2);

      R = cell(1,ncorner);

      RG = speye(2*np);
      opts_rcip =[];

      for icorner=1:ncorner
        clist = corners{icorner}.clist;
        isstart = corners{icorner}.isstart;
        nedge = corners{icorner}.nedge;

        if issymmetric && mod(icorner,2)==0

        else
          rparslocal = [];
          rparslocal.k = k;
          rparslocal.c = c(:,clist);
          rparslocal.coef = coef;

          cparslocal = cell(1,nedge);
          for i=1:nedge
            cparslocal{i} = cpars{clist(i)};
            cparslocal{i}.islocal = isstart(i);
          end

          fcurvelocal = cell(1,nedge);
          for i=1:nedge
            fcurvelocal{i} = @(t) clm.funcurve(t,clist(i),cparslocal{i},icase);
          end

          [Pbc,PWbc,starL,circL,starS,circS,ilist] = rcip.setup(ngl,ndim,nedge,isstart);

          h0 = zeros(1,nedge); % chunk size at the coarsest level
          for i=1:nedge
            if isstart(i)
              h0(i) = chnkr(clist(i)).h(1);
            else
              h0(i) = chnkr(clist(i)).h(end);
            end
          end
          h0 = h0*2; % note the factor of 2 here!!!!

          nsub = 20; % level of dyadic refinement in the forward recursion for computing R
          fkeruse = @(s,ilistl) clm.buildmat_fast(s,rparslocal,opts_rcip,opdims,...
               glwts,ilistl,logquad);
          R{icorner} = rcip.Rcomp_fast_general(ngl,nedge,ndim,Pbc,PWbc,nsub,...
            starL,circL,starS,circS,ilist,...
            h0,isstart,fcurvelocal,glnodes,fkeruse);
        end

        starind = [];
        for i=1:nedge
          if isstart(i)
            starind = [starind inds(clist(i))*ngl*ndim+(1:2*ngl*ndim)];
          else
            starind = [starind inds(clist(i)+1)*ngl*ndim-fliplr(0:2*ngl*ndim-1)];
          end
        end

        M(starind,starind) = 0;

        if issymmetric && mod(icorner,2)==0 
          % due to pointwise block structure, need to reverse the order for
          % each edge, while keeping the same order of 2x2 blocks.
          R11 = R{icorner-1}(1:2:end,1:2:end);
          R12 = R{icorner-1}(1:2:end,2:2:end);
          R21 = R{icorner-1}(2:2:end,1:2:end);
          R22 = R{icorner-1}(2:2:end,2:2:end);

          indinv = [];
          n0 = 2*ngl;
          for i=1:nedge
            indinv = [indinv (i-1)*n0+(n0:-1:1)];
          end

          R11 = R11(indinv,indinv);
          R12 = R12(indinv,indinv);
          R21 = R21(indinv,indinv);
          R22 = R22(indinv,indinv);

          R2 = zeros(2*ngl*nedge*ndim);

          R2(1:2:end,1:2:end) = R11;
          R2(1:2:end,2:2:end) = R12;
          R2(2:2:end,1:2:end) = R21;
          R2(2:2:end,2:2:end) = R22;

          RG(starind,starind) = R2;
          if nonsmoothonly
              M(starind,starind) = inv(R2) - eye(2*ngl*nedge*ndim);
          end
        else
          RG(starind,starind) = R{icorner};
          if nonsmoothonly
              M(starind,starind) = inv(R{icorner}) - eye(2*ngl*nedge*ndim);
          end 
        end
      end

      dt = toc(start);
      disp(['Preconditioner construction time = ', num2str(dt), ' seconds'])
      %fprintf('%5.2f seconds : time to compute R\n',t1)
    end
    
    

end
