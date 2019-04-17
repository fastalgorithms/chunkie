# WARNING

WARNING: chunkie is currently under development and
this readme is somewhat aspirational. Consider it
pre-beta

# chunkie: CHUNK-based Integral Equations

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
- chunkie has routines for setting up system matrices
corresponding to logarithmically singular integral equation
kernels defined on a chunker
- chunkie is designed to inter-operate with Ken Ho's fast
linear algebra in MATLAB package (FLAM) --- more to come!
- chunkie includes various routines for evaluating layer
potentials and functions defined on chunkers

## Installing chunkie

To install chunkie, you need mwrap to compile the
mexfile. With mwrap installed, you should be able
to build the mexfile by executing

```bash
make mexfile
```
in the top-level directory. Once the mexfile is installed,
all that you need to do to use chunkie is to add the mwrap
and matlab directories (with subdirectories) to your matlab
path.

## Using chunkie



## License

chunkie is copyright 2019 Michael O'Neil, James
Bremer, Travis Askham, and Manas Rachh.

chunkie is available under the terms of the
BSD 3-clause license, which should have been included
in the distribution (see LICENSE.md)

chunkie also depends on some routines from Rokhlin,
Greengard, and Gimbutas which were in-house tools that
solve common problems. We have grouped these routines
into the folder external/rgg_tools.

## TO DO

chunkie is new software. Some short term plans
include:

- professionalization of documentation
- building out a suite of demonstrations
- tighter integration with FLAM
- allowing for different singular integration
paradigms (e.g. following the work of Helsing et al.,
Kloeckner et al., and Serkh et al.)
- arbitrary dimension curves on chunks
- libraries for evaluating some common integral
kernels
