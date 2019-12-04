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

Add the chunkie subfolder to your matlab path.
For some features, you will need to also add the
FLAM library to your path. This is included as
a submodule. You can either clone with the submodules

    git clone --recurse-submodules https://github.com/fastalgorithms/chunkie.git

or initialize the submodules after a git pull.
Alternatively, it should work if you already have
a reasonably up-to-date copy of FLAM on your path
see [the FLAM GitHub page](https://github.com/klho/FLAM)


## Using chunkie

Check out demo folder.

## License

chunkie is copyright 2019 the chunkie team

chunkie proper (the contents of the chunkie
folder) is available under the terms of the
BSD 3-clause license, which should have been included
in the distribution (see chunkie/LICENSE.md)

this top-level folder contains some fortran
files and others which may be subject to a slightly
different license.

## chunkie team

chunkers:
- Travis Askham
- Manas Rachh
- Michael O'Neil

kindly donated code:
- Singular quads: James Bremer (chunkie/+chnk/+quad/+brem)
- Classic Fortran routs: Leslie Greengard, Zydrunas
Gimbutas, Vladimir Rokhlin

## TO DO

chunkie is new software. Some short term plans
include:

- professionalization of documentation
- building out a suite of demonstrations
- tighter integration with FLAM
- allowing for different singular integration
paradigms (e.g. following the work of Helsing et al.,
Kloeckner et al., and Serkh et al.)
- libraries for evaluating some common integral
kernels