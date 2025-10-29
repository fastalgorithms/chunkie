
verts = [0 0.3 1 1 .5;
         0 1 1 0 0];
evs = [[1;2],[2;3],[3;4],[4;5],[5;1]];
evs = flipud(evs);
cgrph = chunkgraph(verts,evs);
% figure(1)
% plot(cgrph)
% quiver(cgrph)
% source = [1.5, 1.5]';
% str = []; str.r = source; str.n = [1,0]';
% x = linspace(0,1,120)+1e-2; x = x(1:end-1); y = x;
% [xx,yy] = meshgrid(x,y);
% targs = []; targs.r = [xx(:)';yy(:)'];
% dkern = kernel('lap','d');
% skern = kernel('lap','s');
% spkern = kernel('lap','sp');
% dpkern = kernel('lap','dp');
% 
% kerns(5,5) = kernel();
% for i = 1:3
%     for j = 1:3
%     kerns(i,j) = 2*spkern;
%     end
% end
% for i = 1:3 
%     kerns(4,i) = -2*skern;
%     kerns(i,4) = 2*dpkern;
%     kerns(5,i) = -2*skern;
%     kerns(i,5) = 2*dpkern;
% end
% 
% kerns(4:5,4) = dkern;
% kerns(4:5,5) = dkern;
% 
% 
% chnkg_n = slicegraph(cgrph,1:3);
% chnkg_d = slicegraph(cgrph,4:5);
% A = chunkermat(cgrph,kerns);
% A = A + eye(size(A));
% source = [1.5, 1.5]';
% str = []; str.r = source;
% rhs_n = 2*spkern.eval(str,chnkg_n);
% rhs_d = -2*skern.eval(str,chnkg_d);
% rhs = -[rhs_n; rhs_d];
% b = A\rhs;
% 
% kern_plot(1,5) = kernel();
% for i = 1:3
%     kern_plot(1,i) = skern;
% end
% kern_plot(1,4) = dkern;
% kern_plot(1,5) = dkern;
% 
% 
% int = chunkgraphinregion(cgrph,targs);
% 
% opts = []; opts.eps  = 1e-10;
% sol = chunkerkerneval(cgrph,kern_plot,b,targs,opts);
% u_in = skern.eval(str,targs);
% u_in(int == 1) = NaN;
% u_tot = u_in + sol;
% h = pcolor(xx,yy,reshape(u_tot,size(xx)));
% h.EdgeColor = 'None';
% h = pcolor(xx,yy,reshape(log10(abs(u_tot)),size(xx)));
% h.EdgeColor = 'None';
% colorbar;

%%

isplit = 2;
[cgrph_new, isort] = chnk.fins.split_chunkgraph(cgrph,isplit,3);
figure(3);clf
plot(cgrph_new)
quiver(cgrph_new)
int2 = chunkgraphinregion(cgrph_new,targs);

figure(4);clf
plot(cgrph_new.echnks(isplit))
hold on
plot(cgrph_new.verts(1, cgrph_new.edgesendverts(:,isplit)), cgrph_new.verts(2, cgrph_new.edgesendverts(:,isplit)), 'o')
hold off

[cgrph_new2, isort2] = chnk.fins.split_chunkgraph(cgrph_new,1,3);


norm(cgrph_new.r(:,:) - cgrph.r(:,(isort)))
norm(cgrph_new2.r(:,:) - cgrph_new.r(:,isort2))
norm(cgrph_new2.r(:,:) - cgrph.r(:,isort(isort2)))


%%

iverts = [1,1];
iedges = [1,2];
iedges_ref = 3-iedges;
nchs = [3,2];
cs = [-1, -1];
% iverts = 1;

[cgrph_fin, isort, rfins, idignore] = chnk.fins.build_fins(cgrph, kerns, iverts, iedges, iedges_ref, nchs, cs);
figure(1);clf
plot(cgrph_fin)
% quiver(cgrph_fin)
hold on

norm(cgrph_fin.r(:,:) - cgrph.r(:,isort))

opts =[]; opts.rcip_ignore = [idignore,5];

for i = 1:length(rfins)
    plot(rfins{i}.r(1,:), rfins{i}.r(2,:),'.')
end
axis equal

% 
% 
% 
% [r1, tau1] = chunkends(chnkg.echnks(4),chnkg.echnks(4).nch);
% [r2, tau2] = chunkends(chnkg.echnks(5),1);
% r2 = r2(:,1); tau2 = tau2(:,1);
% r1 = r1(:,2); tau1 = tau1(:,2);
% r1dt1 = dot(r1,tau1);
% r2dt2 = dot(r2,tau2);
% dkern4 = dkern + kernel(@(s,t) dkern.eval(s_flip(s,r1dt1,tau1),t));
% dkern5 = dkern + kernel(@(s,t) dkern.eval(s_flip(s,r2dt2,tau2),t));
% dpkern4 = dpkern + kernel(@(s,t) dpkern.eval(s_flip(s,r1dt1,tau1),t));
% dpkern5 = dpkern + kernel(@(s,t) dpkern.eval(s_flip(s,r2dt2,tau2),t));
% a = s_flip(str,r1dt1,tau1)
% u_in = dkern5.eval(str,targs);
% h = pcolor(xx,yy,reshape(u_in,size(xx)));
% h.EdgeColor = 'None';
% 
% kerns(5,5) = kernel();
% for i = 1:3
%     for j = 1:3
%     kerns(i,j) = 2*spkern;
%     end
% end
% for i = 1:3 
%     kerns(4,i) = -2*skern;
%     kerns(i,4) = 2*dpkern4;
%     kerns(5,i) = -2*skern;
%     kerns(i,5) = 2*dpkern5;
% end
% kerns(1,5).shifted_eval = @(s,t,ori) 2*(dpkern.eval(s,t) + dpkern.eval(s_flip(s,0,tau2),t));
% kerns(3,4).shifted_eval = @(s,t,ori) 2*(dpkern.eval(s,t) + dpkern.eval(s_flip(s,0,tau1),t));
% 
% kerns(4:5,4) = dkern4;
% kerns(4:5,5) = dkern5;
% 
% 
% 
% chnkg_n = slicegraph(chnkg,1:3);
% chnkg_d = slicegraph(chnkg,4:5);
% A = chunkermat(chnkg,kerns);
% A = A + eye(size(A));
% source = [1.5, 1.5]';
% str = []; str.r = source;
% rhs_n = 2*spkern.eval(str,chnkg_n);
% rhs_d = -2*skern.eval(str,chnkg_d);
% rhs = -[rhs_n; rhs_d];
% b = A\rhs;
% 
% kern_plot(1,5) = kernel();
% for i = 1:3
%     kern_plot(1,i) = skern;
% end
% kern_plot(1,4) = dkern4;
% kern_plot(1,5) = dkern5;
% 
% 
% 
% 
% int = chunkgraphinregion(chnkg,targs);
% 
% opts = []; opts.eps  = 1e-10;
% sol = chunkerkerneval(chnkg,kern_plot,b,targs,opts);
% u_in = skern.eval(str,targs);
% u_in(int == 1) = NaN;
% u_tot = u_in + sol;
% h = pcolor(xx,yy,reshape(u_tot,size(xx)));
% h.EdgeColor = 'None';
% h = pcolor(xx,yy,reshape(log10(abs(u_tot)),size(xx)));
% h.EdgeColor = 'None';
% colorbar;
% 
% function src2 = s_flip(src,ndr0,n) 
%     % SCR2 We find the reflected source of src about a line with normal n through r0
%     %
%     %Syntax: scr2 = s_slip(src,ndr0,n)
%     %
%     %Input:
%     %   src - struct holding source location
%     %
%     %   ndr0 - the dot product of the normal vector with the point r_0
%     src2 =[];
%     ndr = n(:)'*src.r(:,:);
%     src2.r = src.r(:,:) - 2*(ndr-ndr0).*n(:);
% 
%     if isfield(src,"n") || isprop(src,"n"), src2.n = src.n; end
% 
% end