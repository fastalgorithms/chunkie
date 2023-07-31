
.. role:: matlab(code)
   :language: matlab   

Solving Standard Boundary Value Problems
=========================================

Chunkie uses the integral equation method to solve PDEs.
In contrast with finite element methods, where the PDE is
typically converted to a weak form, the standard approach for
an integral equation method is to select an integral representation
of the solution that satisfies the PDE inside the domain
*a priori*. This cannot be done for just any PDE; typically,
integral equation methods are limited to constant coefficient,
homogeneous PDEs. But when the methods apply, they are
efficient and accurate.

The following is a very high level discussion of the main ideas.
Some more specific examples are included further below.
For more, see, *inter alia*, [Rok1985]_, [CK2013]_, [GL1996]_, and
[Poz1992]_.

Let :math:`\Omega` be a domain with boundary :math:`\Gamma`.
Suppose that we are interested in the solution of the boundary
value problem

.. math::

   \begin{align*}
   \mathcal{L} u &= 0 & \textrm{ in } \Omega \; ,\\
   \mathcal{B} u &= f & \textrm{ on } \Gamma \; ,
   \end{align*}

where :math:`\mathcal{L}` is a linear, constant coefficient PDE
operator and :math:`\mathcal{B}` is a linear trace operator that
restricts :math:`u`, or some linear operator applied to :math:`u`, to the
boundary.

Suppose that :math:`G` is the Green function for :math:`\mathcal{L}`
and let :math:`K` be an integral kernel that is a linear
operator applied to :math:`G`. Then, :math:`u` can be represented
as a *layer potential*, i.e. it can be written in the form

.. math::

   u(x) = \int_\Gamma K(x,y) \sigma(y) \, ds(y) \; ,

which has the property that :math:`u` automatically solves the PDE
for any choice of density :math:`\sigma` (defined on the boundary
alone). The idea is to then solve for a specific :math:`\sigma`
for which the boundary condition holds:

.. math::

   \mathcal{B} \left [ \int_\Gamma K(\cdot,y)\sigma(y)\, ds(y) \right ]
   = f \textrm{ on } \Gamma \; .

Some of the art in the design of an integral equation method is
in choosing :math:`K` so that the above yields a reasonable boundary
integral equation for :math:`\sigma`. Fortunately, there is a
large literature on good quality integral representations for the most
common problems in classical physics (linear elasticity, electrostatics,
acoustics, and fluid flow). We include some examples below.

This integral equation can then be discretized. Chunkie employs a
Nyström discretization method, i.e. the unknown :math:`\sigma`
is represented by its values at the nodes of a :matlab:`chunker`
object. Because :math:`K` is defined in terms of the Green function
for :math:`\mathcal{L}` it often has mild (or even strong) singularities
when restricted to :math:`\Gamma`. One of the key features of
chunkie is that the function :matlab:`layermat` will automatically return
a Nyström discretization of a layer potential using high-order
accurate quadrature rules [BGR2010]_.

After solving the discrete system for :math:`\sigma`, the PDE solution,
:math:`u`, can then be recovered by evaluating the layer potential
representation at any points of interest. This integral is nearly singular
for points near the boundary. Chunkie provides the function
:matlab:`layereval` for doing this accurately, employing adaptive
integration when necessary.


A Laplace Problem
------------------




A Helmholtz Scattering Problem
-------------------------------



A Stokes Flow Problem
----------------------


References
------------

.. [Rok1985] Rokhlin, Vladimir. "Rapid solution of integral equations
	     of classical potential theory." Journal of Computational
	     Physics 60.2 (1985): 187-207.


.. [Poz1992] Pozrikidis, Constantine. *Boundary integral and singularity
	     methods for linearized viscous flow*. Cambridge University
	     Press, 1992.


.. [GL1996] Guenther, Ronald B., and John W. Lee. *Partial differential
	    equations of mathematical physics and integral equations*.
	    Courier Corporation, 1996.
	    
.. [CK2013] Colton, David, and Rainer Kress. *Integral equation
	    methods in scattering theory*. Society for Industrial
	    and Applied Mathematics, 2013.

.. [BGR2010] Bremer, James, Zydrunas Gimbutas, and Vladimir Rokhlin.
	     "A nonlinear optimization procedure for generalized Gaussian
	     quadratures." SIAM Journal on Scientific Computing 32.4 (2010):
	     1761-1788.
