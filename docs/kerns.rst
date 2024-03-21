.. role:: matlab(code)
   :language: matlab  

The kernel class
==================
The kernel class is an easy way of specifying integral operators
in classical potential theory. 
The kernel object has several attributes including parameters required
to define the kernel, the type of singularity for determining the quadrature
to be used, function handles for evaluating the kernel, 
and fmm acceleration routines if available. 

The constructor for the kernel objects take the form::

    kernel('PDE name', 'kernel name', 'extra params')

Chunkie includes pre-edefined kernels that arise in the solution of the
following PDEs.

- :ref:`lap`
- :ref:`helm`
- :ref:`stokes`
- :ref:`elasticity`

.. _lap:

Laplace kernels
----------------
The free-space Green's function for Laplace's equation in two dimensions,
denoted by $G_{0}(x, y)$, is given by

.. math::

   G_{0}(x,y) = -\frac{1}{2\pi} \log{(\sqrt{(x_{1} - y_{1})^2 + (x_{2} -  y_{2})^2)}} \,,

where $x=(x_{1},x_{2})$, and $y=(y_{1},y_{2})$.

The PDE keywords for using any Laplace kernels are 'Laplace', or 'l'. 
In the following, the variable $y$ will be referred to as the source, and the
variable x will be the target.

The following layer potentials
defined by their associated kernels are supported:

- 'single' or 's': Laplace single layer potential, $G_{0}(x,y)$
- 'double' or 'd': Laplace double layer potential, $n(y) \cdot \partial_{y} G_{0}(x,y)$
- 'combined' or 'c': Laplace combined layer potential, $c_{1} n(y) \cdot
  \partial_{y} G_{0}(x,y) + c_{2} G_{0}(x,y)$, where the coefficients
  $c_{1},c_{2}$ are passed as an extra array of length 2, with default values of
  $(c_{1},c_{2}) = (1,1)$.
- 'sprime' or 'sp': Normal derivative of Laplace single layer
   potential, $n(x) \cdot \partial_{x} G_{0}(x,y)$
- 'stau' or 'st': Tangential derivative of Laplace single layer potential,
   $\tau(x) \cdot \partial_{x} G_{0}(x,y)$
- 'dprime' or 'dp': Normal derivative of Laplace double layer potential, 
   $n(x) \cdot \partial_{x} n(y) \cdot \partial_{y} G_{0}(x,y)$
- 'sgrad' or 'sg': Gradient of Laplace single layer potential, 
   $\partial_{x} G_{0}(x,y)$
- 'dgrad' or 'dg': Gradient of Laplace double layer potential,
   $\partial_{x} n_{y} \cdot \partial_{y} G_{0}(x,y)$
