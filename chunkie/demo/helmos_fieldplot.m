function helmos_fieldplot(ufield, xg, yg, varargin)
%HELMOS_FIELDPLOT  Plot a field on a regular tensor grid (testhelmos style).
%
%   helmos_fieldplot(ufield, xg, yg) plots abs(ufield) reshaped over xg,yg.
%
%   helmos_fieldplot(ufield, xg, yg, 'mode', m) where m is one of:
%       'abs'   (default)  imagesc of |u|
%       'real'              imagesc of real(u) symmetric scale
%       'log10err'          imagesc of log10|u|, useful when ufield is the
%                           difference between numerical and reference
%
%   Optional name-value pairs:
%       'cg'        chunkgraph to overlay (plotted as black curves)
%       'xylim'     [xmin xmax ymin ymax] axis limits (default tight to xg,yg)
%       'title'     figure title
%
%   Convention: ufield is a column vector of length ngr*ngr where ufield(k)
%   corresponds to the target with linear index k = (j-1)*ngr + i, i.e.
%   columns of (xg,yg) are stacked. ngr = length(xg) = length(yg).

p = inputParser;
p.addParameter('mode','abs');
p.addParameter('cg',[]);
p.addParameter('xylim',[]);
p.addParameter('title','');
p.parse(varargin{:});
mode = p.Results.mode;
cg   = p.Results.cg;
xylim = p.Results.xylim;
ttl   = p.Results.title;

ngr = length(xg);
assert(length(yg) == ngr, 'helmos_fieldplot: xg and yg must have same length');
assert(numel(ufield) == ngr*ngr, 'helmos_fieldplot: ufield length must be ngr*ngr');

F = reshape(ufield(:), ngr, ngr);

figure;
set(gca,'FontSize',12);
hold on;
switch lower(mode)
    case 'abs'
        imagesc(xg, yg, abs(F));
        colormap(jet);
    case 'real'
        m = max(abs(F(:)));
        imagesc(xg, yg, real(F), [-m m]);
        colormap(jet);
    case 'log10err'
        L = log10(abs(F));
        m = max(L(:));
        imagesc(xg, yg, L, [log10(eps) m]);
        colormap(flipud(pink(256)));
    otherwise
        error('helmos_fieldplot: unknown mode "%s"', mode);
end
axis xy;
colorbar;
xlabel('$x$','Interpreter','LaTeX','FontSize',17);
ylabel('$y$','Interpreter','LaTeX','FontSize',17);

if ~isempty(cg)
    for ie = 1:length(cg.echnks)
        e = cg.echnks(ie);
        rr = reshape(e.r, 2, []);
        plot(rr(1,:), rr(2,:), 'k-', 'LineWidth', 1.0);
    end
end

if ~isempty(ttl), title(ttl); end
axis equal;
if isempty(xylim)
    xylim = [min(xg) max(xg) min(yg) max(yg)];
end
axis(xylim);
drawnow;
end
