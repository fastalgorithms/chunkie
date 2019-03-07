# chunkie --- CHUNK-based Integral Equations

A MATLAB package for prototyping integral equation
methods in two dimensions.

While chunkie is primarily intended as a proto-typing
or pedagogical tool, it is designed to be reasonably
efficient and has been used to produce research-grade
results.

At-a-glance:
- given a parametrization of a curve, chunkie will return
a "chunker" object which stores the description of the
curve in chunk format (the curve is discretized into chunks
such that on each chunk a legendre interpolant in parameter
space is accurate to some prescribed accuracy).
- routines for setting up system matrices corresponding
to logarithmically singular integral equation kernels
on a chunker
- designed to inter-operate with Ken Ho's fast linear algebra
in MATLAB (FLAM, more to come!)
- various routines for evaluating layer potentials
and functions defined on chunkers

## 

