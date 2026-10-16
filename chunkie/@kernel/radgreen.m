function obj = radgreen(green, type, varargin)
%KERNEL.RADGREEN   Construct a kernel from a radial Green's function.
%   KERNEL.RADGREEN(GREEN, TYPE) constructs the kernel of the given TYPE
%   for the radial Green's function G(|x-y|) given by the function handle
%   GREEN, which is expected to be in the format of CHNK.LAP2D.GREEN, with
%   derivatives with respect to the target.
%
%   KERNEL.RADGREEN(GREEN, TYPE, COEFS) or KERNEL.RADGREEN(GREEN, TYPE, NU)
%   sets the parameter of the types which take one.
%
%   KERNEL.RADGREEN(GREEN, TYPE, ..., 'dim', DIM) sets the dimension of
%   the points and derivative tensors, 2 (default) or 3. The normals,
%   tangents and curvatures of the plate kernels are only defined for 2.
%
%   KERNEL.RADGREEN(GREEN, TYPE, ..., 'sing', SING) sets the singularity
%   type of the kernel, 'log', 'pv' or 'hs' (see KERNEL). By default it is
%   found from the decay of the kernel on a test curve.
%
%   The types, with D = G_{n_y}, S = G, D' = G_{n_x n_y}, S' = G_{n_x}:
%
%      TYPE                                       KERNEL
%      ----                                       ------
%      's' ('single')                             S
%      'd' ('double')                             D
%      'sp' ('sprime')                            S'
%      'dp' ('dprime')                            D'
%      'c' ('combined')                           COEFS(1)*D + COEFS(2)*S
%      'cp' ('cprime')                            COEFS(1)*D' + COEFS(2)*S'
%      'sg' ('sgrad')                             grad S, opdims [dim 1]
%      'dg' ('dgrad')                             grad D, opdims [dim 1]
%      'cg' ('cgrad')                             COEFS(1)*grad D +
%                                                 COEFS(2)*grad S
%      'st' ('stau')                              tangential derivative of S
%      'hilb' ('hilbert')                         two times the adjoint of
%                                                 'st', the Hilbert
%                                                 transform for Laplace
%      's2trans'                                  [S; S']
%      'd2trans'                                  [D; D']
%      'c2trans' ('c2t')                          [C; C'] with C the
%                                                 combined layer, COEFS is
%                                                 2x2, as in KERNEL.HELM2D
%      'trans' ('all', 'trans_sys', 'tsys')       [COEFS(1,1)*D COEFS(1,2)*S;
%                                                  COEFS(2,1)*D' COEFS(2,2)*S']
%      'trans_rep' ('trep')                       [COEFS(1)*D COEFS(2)*S]
%      'trans_rep_prime' ('trep_p', 'trans_rep_p') [COEFS(1)*D' COEFS(2)*S']
%      'trans_rep_grad' ('trep_g', 'trans_rep_g') [COEFS(1)*grad D
%                                                  COEFS(2)*grad S]
%
%   and the plate kernels, for dim = 2, with Poisson ratio NU:
%
%      'clamped_bc' ('clamped_plate_bcs')         boundary conditions of
%                                                 the clamped plate on a
%                                                 point source, [G; G_{n_x}]
%      'clamped_eval' ('clamped_plate_eval')      kernels of the clamped
%                                                 plate representation
%      'supported_bc' ('supported_plate_bcs')     boundary conditions of
%                                                 the supported plate, NU
%      'supported_eval' ('supported_plate_eval') kernels of the supported
%                                                 plate representation, NU,
%                                                 srcinfo.data(1,:) is the
%                                                 derivative of curvature
%      'free_bc' ('free_plate_bcs')               boundary conditions of
%                                                 the free plate, NU
%      'free_eval' ('free_plate_eval')            kernels of the free plate
%                                                 representation, NU
%
%   The plate kernels are those of CHNK.FLEX2D.KERN with the Green's
%   function replaced by GREEN. The integral equation kernels of
%   CHNK.FLEX2D.KERN (clamped_plate, free_plate, supported_plate) have
%   terms specific to the flexural wave Green's function and are not
%   included. See CHNK.RADPP.KERN for the definitions.
%
%   KERNEL.RADGREEN(GREEN, FUN, 'nderiv', ND, 'opdims', [M N]) constructs
%   a custom kernel. FUN is a function handle, W = FUN(D, SRCINFO,
%   TARGINFO), which returns the (nt,ns,M*N) entries of the kernel block,
%   in column-major order, as combinations of derivatives of G. D(TD, SD)
%   is the derivative of G along the vectors in the cell array TD, with
%   respect to the target, and those in SD, with respect to the source.
%   TD vectors are (dim,nt) arrays and SD vectors are (dim,ns) arrays (or
%   (dim,1) for a constant direction). At most ND <= 4 derivatives in
%   total are used. For example, the combined layer kernel is
%
%      fun = @(D,s,t) coefs(1)*D({},{s.n}) + coefs(2)*D({},{});
%      K = kernel.radgreen(green, fun, 'nderiv', 1);
%
%   Example:
%      zk = 1.1;
%      green = @(s,t) chnk.helm2d.green(zk, s, t);
%      K = kernel.radgreen(green, 'c', [1 1i]);
%
% See also KERNEL.RADPP, CHNK.RADPP.KERN, CHNK.RADPP.SUMGREENS.

if ( nargin < 1 )
    error('Missing Green''s function.');
end
if ( nargin < 2 )
    type = 's';
end

obj = chnk.radpp.build(green, [], type, varargin, false);

end
