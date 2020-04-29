.. chunkie documentation master file, created by
   sphinx-quickstart on Tue Apr 28 21:16:39 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


chunkie
========

.. toctree::
   :maxdepth: 2
   :caption: Contents:


chunkie is a MATLAB package for solving partial differential
equations in two dimensional domains. The name "chunkie" derives
from the use of chunks (panels with scaled Legendre nodes) to
discretize the domain boundary and the use of an integral equation
method to solve the PDE.

With chunkie, you can solve the Helmholtz equation on an exterior
domain in a few lines of MATLAB:

    import project
    # Get your stuff done
    project.do_stuff()

Features
--------

- easy-to-use: includes pre-packaged solvers for Laplace, Helmholtz,
  and Stokes equations and convenient tools for visualizing
  domains and plotting solutions
- flexible: allows you to define your own integral kernels and
  use the quadrature routines in a modular fashion
- fast: works with FLAM (fast linear algebra in MATLAB)
  for fast, direct integral equation solves

Installation
------------

chunkie is pure MATLAB. You can download ...

Contributing
------------

We welcome your collaboration! First, see the
`developer guide <https://github.com/fastalgorithms/chunkie/wiki>`_.

- `Issue Tracker <https://github.com/fastalgorithms/chunkie/issues>`_
- `Source Code <https://github.com/fastalgorithms/chunkie>`_


License
-------

The project is licensed under a modified BSD 3-clause
license.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
