function chnkr = roundedpolygon(verts,cparams,pref,edgevals)
%CHUNKPOLY return a chunker object corresponding to the polygon with the
% given vertices. By default, a polygon with rounded corners is returned.
% Open and higher dimensional "polygons" are allowed. Optional dyadic
% refinement near corners of true polygons (no rounding).
%
% Syntax: chnkr = chunkpoly(verts,cparams,pref,edgevals)
%
% Input:
%    verts - (dimv,nverts) array of "polygon" vertices
%            in order
%
% Optional input:
%    cparams - options structure
%       cparams.rounded = true if corner rounding is to be used.
%                         false if no rounding is used (true)
%       cparams.widths = radius around each corner point where either the
%                        rounding or dyadic refinement is applied 
%                        (defaults to 1/10th of minimum
%                         of length of adjoining edges)
%       cparams.autowidths = automatically compute widths (false)
%       cparams.autowidthsfac = if using autowidths, set widths
%                             to autowidthsfac*minimum of adjoining
%                             edges (0.1)
%    	cparams.ifclosed = true, if it's a closed polygon
%                          false, if it's an open segment (true)
%       cparams may also specify any of the values to send to the
%            refine routine, otherwise it is called with defaults
%
%           ~ Rounding parameters ~
% 	    cparams.eps - resolve curve to tolerance eps
%                    resolve coordinates, arclength,
%          	     and first and second derivs of coordinates
%		         to this tolerance (1.0e-6)
%
%   pref - chunkerpref object/ preference structure 
%       pref.k determines order of underlying Legendre nodes
%   edgevals - very optional input. specifies constant values along each
%              edge. The routine then stores a smoothed interpolant 
%              of these edge values on the rounded structure in the 
%              output chunker's data field
%
% Output:
%   chnkr - chunker object corresponding to (rounded) polygon
%                  
% Examples:
%   barbell_verts = chnk.demo.barbell();
%   chnkr = roundedpolygon(barbell_verts); % rounded "barbell" domain
%                                     % with standard options
%   cparams = []; cparams.rounded = false;
%   pref = []; pref.k = 30;
%   chnkr = chunkpoly(barbell_verts,cparams,pref); % not rounded
%
% See also CHUNKERFUNC, CHUNKER, CHUNKERPREF, REFINE

chnkr = chunkerpoly(verts,cparams,pref,edgevals);