

.. role:: matlab(code)
   :language: matlab   

Describing More Complicated Geometries with Chunkgraphs
========================================================

Many practical applications, especially those involving corners and 
multiple junction interfaces, require a more detailed description of the 
geometry, beyond an array of `chunker` objects. While it is relatively
easy to specify these more complicated problems, integral equation methods
require special care in the vicinity of corners 
and mutltiple junctions in order to achieve high precision. 
The requisite additional information can be
stored in `chunkie` in a `chunkgraph` object, which is essentially
a collection of `chunker` objects and a graph-theoretic description
of their connectivity. In particular, each smooth boundary component,
stored as a `chunker`, corresponds to an edge of the graph, and
points where multiple smooth boundary components meet (or where a
smooth component ends), correspond to vertices of the graph.

Creating Chunkgraphs
--------------------
The easiest way to define :matlab:`chunkgraphs` is using its constructor
which requires at least two inputs, a collection of vertices and
connectivity information specified as a :matlab:`(2,nedges)` array
with each column corresponding to the indices of the vertices
defining the edge. The edge is directed and begins at the vertex
in the first row and ends at the vertex in the corresponding second row.

In the example below, we construct a polygon defined as a chunkgraph. As 
before we can visualize the :matlab:`chunkgraph` using overloaded MATLAB
commands :matlab:`plot` and :matlab:`quiver`.

.. include:: ../../chunkie/guide/guide_chunkgraphs.m
   :literal:
   :code: matlab
   :start-after: % START POLYGON
   :end-before: % END POLYGON
	  
.. image:: ../../chunkie/guide/guide_chunkgraphs_polygon.png
   :width: 350px
   :alt: polygon chunkgraph 
   :align: center


Instead of having straight edges, the user can define curved edges
between vertices as well. This is specified either as a single function
handle, in which case all edges would be defined by the same function handle,
or via a cell array of function handles, where vertices are connected
by a straight edge if the corresponding element of the cell 
array of function handles is empty. Each function handle
must have the same specification as the ones used to define :matlab:`chunkers`,
i.e., if :matlab:`fcurve` is a MATLAB function defining the curve
parametrization and :matlab:`t` is a vector fo points in parameter space,
then the output of :matlab:`fcurve` is of the form :matlab:`[r] = fcurve(t)`
or :matlab:`[r,d] = fcurve(t)` or :matlab:`[r,d,d2] = fcurve(t)` where
:matlab:`r` are the coordinates of the points on the curve, :matlab:`d` is
the first derivative of the position with respect to the underlying 
parametrization and :matlab:`d2` is the corresponding second derivative.
For ease of specification, the curves described by the function handles
aren't required to terminate exactly at the appropriate vertices;
the :matlab:`chunkgraph` constructor
will snap the edge to its vertices via an appropriate affine transformation.

Here is an example illustrating a square whose edges are defined via a sine
wave.

.. include:: ../../chunkie/guide/guide_chunkgraphs.m
   :literal:
   :code: matlab
   :start-after: % START CURVED POLYGON
   :end-before: % END CURVED POLYGON
	  
.. image:: ../../chunkie/guide/guide_chunkgraphs_curved_edges.png
   :width: 350px
   :alt: curved polygon chunkgraph 
   :align: center



.. note::
   
   The default domain of parameter space defining 
   an edge of a :matlab:`chunkgraph` is [0,1] as opposed
   to that of a :matlab:`chunker` which is [0,2\pi).
   
   
The same constructor interface is capable of creating complicated 
domains including ones with loops starting and ending at the same vertex,
multiple junctions, and nested regions. The constructor is also capable
of handling :matlab:`chunkgraphs` with smooth closed curves which are
treated as edges with no vertices. All of these are illustrated 
in the construction of the complicated ``hawaiian earing`` with a nested
star-fish below.

.. include:: ../../chunkie/guide/guide_chunkgraphs.m
   :literal:
   :code: matlab
   :start-after: % START HAWAIIAN EARING 
   :end-before: % END HAWAIIAN EARING
	  
.. image:: ../../chunkie/guide/guide_chunkgraphs_hawaiian.png
   :width: 350px
   :alt: Hawaiian earing with inclusion chunkgraph 
   :align: center


Working with chunkgraphs
------------------------
Most functions overloaded for chunkers also work with chunkgraphs. 
As with :matlab:`chunkers`,users are free to edit the position and derivative
fields of the :matlab:`chunkgraph` objects, or even replace the individual
:matlab:`chunkers` that define the :matlab:`chunkgraph`, and the software will not
check if the user has updated the :matlab:`chunkgraph` in a consistent manner. 

In addition to these routines, we provide an additional plotting routine 
which labels the edges and regions defined by the :matlab:`chunkgraph`. 
This functionality is illustrated in the example below.

.. include:: ../../chunkie/guide/guide_chunkgraphs.m
   :literal:
   :code: matlab
   :start-after: % START REGION TEST 
   :end-before: % END REGION TEST 
	  
.. image:: ../../chunkie/guide/guide_chunkgraphs_regions.png
   :width: 350px
   :alt: chunkgraph region plot 
   :align: center

For more details of how this information is stored in the 
:matlab:`chunkgraph` object and to see a survey of other available 
methods, see the :matlab:`chunkgraph` class 
documentation:

.. include:: ../../chunkie/@chunkgraph/chunkgraph.m
   :literal:
   :code: matlab
   :start-after: classdef chunkgraph 
   :end-before: % Syntax:

.. note::

   To obtain the documentation of a class method which has
   overloaded the name of a MATLAB built-in, use the syntax
   :matlab:`help class_name/method_name`. For example:

   .. code:: matlab

      help chunkgraph/move
