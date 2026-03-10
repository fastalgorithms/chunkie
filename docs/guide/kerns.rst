.. role:: matlab(code)
   :language: matlab  

The Kernel Class
==================
The :matlab:`kernel` class provides a convenient way to specify integral
operators in classical potential theory. 
The kernel object has several attributes including parameters required
to define the kernel, the type of singularity for determining the quadrature
to be used, function handles for evaluating the kernel, 
and fast multipole method (FMM) acceleration routines if available. 

The constructor for the kernel objects take the form

.. code:: matlab

    kernel('PDE name', 'kernel name', varargin)


The kernel class documentation lists the information that can be stored
in the object
    
.. literalinclude:: ../../chunkie/@kernel/kernel.m
   :language: matlab
   :start-after: classdef kernel
   :end-before: % author

The :matlab:`eval` property and the :matlab:`ptinfo` format
---------------------------------------------------------------

The :matlab:`eval` property in the kernel class is the function which
actually evaluates the interaction between sources and targets. It
expects two structs, usually denoted :matlab:`s` and :matlab:`t`,
which specify the sources and targets respectively.

Points are specified as structs, with the same naming conventions
as the properties of a :matlab:`chunker` object. We refer to this
specification as the :matlab:`ptinfo` format. In particular, if
:matlab:`pts` is a collection of points in the :matlab:`ptinfo` format,
then it has the following properties:

- pts.r - 2 x npts array of (x,y) coordinates for each point [required]
- pts.d - 2 x npts array of (x,y) coordinates of derivative for each point [optional]
- pts.d2 - 2 x npts array of (x,y) coordinates of second derivative for each point [optional]
- pts.n - 2 x npts array of (x,y) normal vector for each point [optional]
- pts.data - m x npts array of m-vectors of data for each point [optional, not fully supported]

Aside from the :matlab:`pts.r` property, the remaining properties are optional
and most are only meaningful for points on a curve.

Suppose that :matlab:`kern` is the kernel object corresponding to a
scalar kernel function of the form $K(x,y)$ and that :matlab:`s` and
:matlab:`t` are structs specifying a collection of :matlab:`ns` sources and
:matlab:`nt` targets, respectively. Then, :matlab:`kern.eval(s,t)` will
return an :matlab:`nt x ns` matrix with the entries

.. math::

   \begin{bmatrix}
   K(x^{(1)},y^{(1)}) & K(x^{(1)},y^{(2)}) & \cdots & K(x^{(1)},y^{(\texttt{ns})}) \\
   K(x^{(2)},y^{(1)}) & K(x^{(2)},y^{(2)}) & \cdots & K(x^{(2)},y^{(\texttt{ns})}) \\
   \vdots & \vdots & & \vdots \\
   K(x^{(\texttt{nt})},y^{(1)}) & K(x^{(\texttt{nt})},y^{(2)}) & \cdots & K(x^{(\texttt{nt})},y^{(\texttt{ns})}) 
   \end{bmatrix}

where $x^{(i)}$ is the $i$th target with coordinates :matlab:`t.r(:,i)` and
$y^{(j)}$ is the $j$th source with coordinates :matlab:`s.r(:,j)`.

If $K$ is a vector- or matrix-valued kernel, then the dimensions of the kernel
are stored in the :matlab:`opdims` property and the output is a
:matlab:`(nt*kern.opdims(1)) x (ns*kern.opdims(2))` matrix. The matrix
is organized into :matlab:`nt x ns` blocks, each of size
:matlab:`kern.opdims(1) x kern.opdims(2)`. For example, if $K$ is a $2 \times 2$
matrix valued kernel, then :matlab:`opdims = [2 2]` and :matlab:`kern.eval(s,t)` will
return a :matlab:`(nt*kern.opdims(1)) x (ns*kern.opdims(2))` matrix with the
entries

.. math::

   \begin{bmatrix}
   K_{11}(x^{(1)},y^{(1)}) & K_{12}(x^{(1)},y^{(1)}) & \cdots & K_{11}(x^{(1)},y^{(\texttt{ns})}) & K_{12}(x^{(1)},y^{(\texttt{ns})}) \\
   K_{21}(x^{(1)},y^{(1)}) & K_{22}(x^{(1)},y^{(1)}) & \cdots & K_{21}(x^{(1)},y^{(\texttt{ns})}) & K_{22}(x^{(1)},y^{(\texttt{ns})}) \\
   \vdots & \vdots & & \vdots & \vdots \\
   K_{11}(x^{(\texttt{nt})},y^{(1)}) & K_{12}(x^{(\texttt{nt})},y^{(1)}) & \cdots & K_{11}(x^{(\texttt{nt})},y^{(\texttt{ns})}) & K_{12}(x^{(\texttt{nt})},y^{(\texttt{ns})}) \\
   K_{21}(x^{(\texttt{nt})},y^{(1)}) & K_{22}(x^{(\texttt{nt})},y^{(1)}) & \cdots & K_{21}(x^{(\texttt{nt})},y^{(\texttt{ns})}) & K_{22}(x^{(\texttt{nt})},y^{(\texttt{ns})}) \\
   \end{bmatrix}

  
Defining your own kernel class object
--------------------------------------

While the integral kernels for many common PDEs are available
as built-ins in chunkIE, the utilities will work for a broad range
of integral kernels. 
To define your own kernel class object you can first obtain an
empty kernel via :matlab:`kern=kernel()` and then update the
properties as needed. It is generally simpler to define the :matlab:`eval`
function in a separate file. The example below implements the kernel
evaluator for the kernel $K(x,y) = \exp(\imath k|x-y|)$ in a file called
``expkernel.m``. 

.. literalinclude:: ../../chunkie/guide/expkernel.m
   :language: matlab

You can then create an empty kernel object and point the eval
property to this function:

.. literalinclude:: ../../chunkie/guide/guide_kernels.m
   :language: matlab
   :start-after: % START CREATE KERNEL
   :end-before: % END CREATE KERNEL		


Combining kernel objects
--------------------------

Given two kernel objects with compatible dimensions, it is
possible to add and scale them.

.. literalinclude:: ../../chunkie/guide/guide_kernels.m
   :language: matlab
   :start-after: % START ADDING KERNELS
   :end-before: % END ADDING KERNELS

For a :matlab:`chunkgraph` with :matlab:`ne` edges, it is
possible to specify interactions between each pair of edges
in a :matlab:`ne x ne` matrix of kernel objects. Note that
most chunkIE routines will require that
the :matlab:`opdims(2)` property is the same for
all kernels within any column of the kernel matrix and 
the :matlab:`opdims(1)` property is the same for
all kernels within any row of the kernel matrix.

A simple way to create a matrix of kernel objects is to
first build an empty matrix of the desired size and then
to assign each entry a kernel.


.. literalinclude:: ../../chunkie/guide/guide_kernels.m
   :language: matlab
   :start-after: % START MATRIX OF KERNELS
   :end-before: % END MATRIX OF KERNELS


A matrix of kernels object can also be used to
create matrix-valued kernels. A matrix of kernels can
be merged into a single matrix-valued kernel using the
:matlab:`kernel` constructor.

.. literalinclude:: ../../chunkie/guide/guide_kernels.m
   :language: matlab
   :start-after: % START INTERLEAVE
   :end-before: % END INTERLEAVE

	
Built-in Kernels  
---------------------

chunkIE literalincludes built-in kernels that arise in the solution of the
following PDEs.

- :ref:`lap`
- :ref:`helm`
- :ref:`stokes`
- :ref:`elasticity`

Some common notation used below is that $x=(x_1,x_2)$ will refer to
a "target" location and $y=(y_1,y_2)$ will refer to a "source location.
Likewise, $n(x)$ refers to the curve normal at the target and $n(y)$ refers
to the normal at the source. The notation $|x-y|$ denotes the distance between
the source and target, i.e.

.. math::

   |x-y| = \sqrt{ (x_1-y_1)^2 + (x_2-y_2)^2 } \; .


.. _lap:

Laplace kernels
~~~~~~~~~~~~~~~~

The free-space Green's function for Laplace's equation in two dimensions,
denoted by $G_{0}(x, y)$, is given by

.. math::

   G_{0}(x,y) = -\frac{1}{2\pi} \log |x-y| \,,

where $x=(x_{1},x_{2})$, and $y=(y_{1},y_{2})$.

The kernel $S(x,y)=G_0(x,y)$
is the integral kernel for the single layer potential. The double layer kernel is
then $D(x,y) = n(y)\cdot \nabla_y G_0(x,y)$. All Laplace integral kernels are
derived from these two.

For each of the single layer, double layer,
and combined layer (a linear combination of single and double layer), there
are also kernels for the normal derivative (with the name "prime"), tangential
derivative (with the name "tau"), and the gradient (with the name "grad") taken
at the target location. For a kernel $K(x,y)$, these are the kernels
$n(x)\cdot \nabla_xK(x,y)$, $\tau(x)\cdot \nabla_x K(x,y)$, and
$\nabla_x K(x,y)$, respectively.

The PDE name for using Laplace kernels can be 'Laplace' or 'l'. These kernel
types can be obtained with either calling sequence below

.. code:: matlab
	  
	  kern = kernel('l',type)
	  kern = kernel.lap2d(type)

The documentation of the :matlab:`kernel.lap2d` function has details on any
expected additional arguments and default values:

.. literalinclude:: ../../chunkie/@kernel/lap2d.m
   :language: matlab
   :end-before: % author

.. _helm:

Helmholtz kernels
~~~~~~~~~~~~~~~~~~~

The free-space Green's function for the Helmholtz equation in two dimensions,
denoted by $G_{k}(x, y)$, is given by

.. math::

   G_{k}(x,y) = \frac{\imath}{4} H_0^{(1)}(k |x-y|) \, ,

where $k$ is the wavenumber, $x=(x_{1},x_{2})$, and $y=(y_{1},y_{2})$. 

The documentation of the :matlab:`kernel.helm2d` function has details on any
expected additional arguments and default values:

.. literalinclude:: ../../chunkie/@kernel/helm2d.m
   :language: matlab
   :end-before: % author


.. _stokes:

Stokes kernels
~~~~~~~~~~~~~~~~

The free-space Green's function for the Stokes equation in two dimensions,
known as a Stokeslet and denoted by ${ G}(x, y)$, is a $2\times 2$
tensor whose $i,j$ entry is

.. math::

   G_{ij}(x,y) = \frac{1}{4\pi \mu} \left ( \frac{(x_i-y_i)(x_j-y_j)}{|x-y|^2}
   - \delta_{ij} \log |x-y| \right ) \, ,

where $\mu$ is the dynamic viscosity, $x=(x_{1},x_{2})$, and $y=(y_{1},y_{2})$.
The quantity ${ G} \cdot { f}$ gives the velocity induced
by a point charge at $y$ with strength and orientation determined by
$f = (f_1,f_2)$. The corresponding pressure is ${P}\cdot { f}$
where

.. math::

   P_{i}(x,y) = \frac{x_i-y_i}{2\pi |x-y|^2} \, .

Let ${ \sigma}({ u},p)$ denote the Cauchy stress tensor for 
given velocity and pressure fields, i.e.

.. math::

   { \sigma}({ u},p) = 2 \mu { \epsilon}
   + p I = \mu ( \nabla { u} + \nabla { u}^T) + p I \; ,

where $\epsilon$ is the strain rate.

On a curve in the fluid, the corresponding traction vector is
$t(x) = n(x) \cdot \sigma(u,p)$. Specifying the traction is
a common boundary condition for Stokes problems. Kernels corresponding
to the traction of some base kernel are denoted by the suffix
"trac", e.g. :matlab:`strac` is the name for the traction of the
Stokeslet. For reference, the 3-tensor corresponding to the
stress of each column of the Stokeslet, denoted by $T$ and
known as a "stresslet", has the entries

.. math::

   T_{ijk}(x,y) = -\frac{(x_i-y_i)(x_j-y_j)(x_k-y_k)}{\pi |x-y|^4} \; ,

so that the traction of a Stokeslet is $n(x) \cdot T(x,y) \cdot f(y)$.   

The Stokes double layer kernel is defined in terms of the stresslet
as follows

.. math::

   D_{ij}(x,y) = -T_{jki}(x,y) n_k(y) \; .


.. literalinclude:: ../../chunkie/@kernel/stok2d.m
   :language: matlab
   :end-before: % author

.. _elasticity:

Linear elasticity kernels
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :ltrim:

The Lam |eacute| parameters $\lambda$ and $\mu$ are used to specify
the stress-strain relations in linear elasticity, i.e.

.. math::

   \sigma = \lambda \operatorname{Trace} (\epsilon) I + 2\mu \epsilon \; .

The governing equations for a linear elastic solid, in the absence of
body forces, are $\nabla \cdot \sigma = 0$. The displacement induced by
a point force with strength and orientation given by $f$ is
$G \cdot f$, where

.. math::

   G_{ij}(x,y) = \frac{1}{4\pi \mu(\lambda + 2\mu)} \left (
   (\lambda + 3\mu) \delta_{ij} \log|x-y| - (\lambda + \mu)
   \frac{(x_i-y_i)(x_j-y_j)}{|x-y|^2} \right ) \; .

The standard traction and double layer operators are analogous to the
Stokes case.


.. literalinclude:: ../../chunkie/@kernel/elast2d.m
   :language: matlab
   :end-before: % author

   
