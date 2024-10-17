
Describing More Complicated Geometries with Chunkgraphs
========================================================

In many practical applications, espcially those involving corners and 
multiple junction interfaces, require a more complicated description of the 
geometry beyond `chunkers`. Apart from the ease of specification of such complicated
problems, integral equation methods require special care in the vicinity of corners 
and mutltiple junctions in order to achieve high precision. 
Such complicated problems which store the requisite additional information can be
stored in a `chunkgraph` object, which is essentially a collection of `chunkers`
connected between a set of vertices.

Polygon
Polygon with curved edges
Multiple junction interfaces


Working with chunkgraph
Most functions overloaded for chunkers also work with chunkgraph.

Additional region handling capabilities


