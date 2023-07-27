
Curve Discretization with Chunkers
===================================

In chunkie, a smooth, regular curve is discretized by diving
it into pieces, called "chunks", which are then represented
by a polynomial interpolant at scaled Legendre nodes. This
information is stored in a chunker object. 

A chunker object can be obtained from a curve parameterization
using the function chunkerfunc. The chunker class overloads
some MATLAB commands, like plot and quiver, to simplify common
visualization tasks. The code below creates a chunker object
for a circle and plots the geometry and the normal vectors.

.. code:: matlab
	  
   % chunk up circle
	  
   r = 0.7; ctr = [1.0;-0.5];
   circfun = @(t) ctr + r*[cos(t(:).');sin(t(:).')];
   
   chnkr = chunkerfunc(circfun); 

   % plot curve, nodes, and normals
   
   figure()
   plot(chnkr,'r-x')
   hold on
   quiver(chnkr,'b')

.. image:: ../assets/images/guide/chunkers_circle.png
   :width: 350px
   :alt: circle chunker
   :align: center



