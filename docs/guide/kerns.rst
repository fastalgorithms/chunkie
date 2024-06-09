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
    
.. include:: ../../chunkie/@kernel/kernel.m
   :literal:
   :code: matlab
   :start-after: classdef kernel
   :end-before: % author

The :matlab:`ptinfo` format
----------------------------



Built-in Kernels  
---------------------

Chunkie includes built-in kernels that arise in the solution of the
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
then $D(x,y) = n(y)\\cdot \\nabla_y G_0(x,y)$. All Laplace integral kernels are
derived from these two.

For each of the single layer, double layer,
and combined layer (a linear combination of single and double layer), there
are also kernels for the normal derivative (with the name "prime"), tangential
derivative (with the name "tau"), and the gradient (with the name "grad") taken
at the target location. For a kernel $K(x,y)$, these are the kernels
$n(x)\\cdot \\nabla_xK(x,y)$, $\\tau(x)\\cdot \\nabla_x K(x,y)$, and
$\\nabla_x K(x,y)$, respectively.

The PDE name for using Laplace kernels can be 'Laplace' or 'l'. These kernel
types can be obtained with either calling sequence below

.. code:: matlab
	  
	  kern = kernel('l',type)
	  kern = kernel.lap2d(type)

The documentation of the :matlab:`kernel.lap2d` function has details on any
expected additional arguments and default values:

.. include:: ../../chunkie/@kernel/lap2d.m
   :literal:
   :code: matlab
   :start-line: 1
   :end-before: % author

.. _helm:

Helmholtz kernels
~~~~~~~~~~~~~~~~~~~

.. _stokes:

Stokes kernels
~~~~~~~~~~~~~~~~

.. _elasticity:

Linear elasticity kernels
~~~~~~~~~~~~~~~~~~~~~~~~~~


Defining your own kernel class object
--------------------------------------



Combining kernel objects
--------------------------


Defining a matrix of kernel objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Adding and scaling kernel objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


   
