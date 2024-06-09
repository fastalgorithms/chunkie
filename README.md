# chunkIE: a MATLAB integral equation toolbox

A MATLAB package for prototyping integral equation
methods in two dimensions.
While chunkie is primarily intended as a proto-typing
or pedagogical tool, it is designed to be reasonably
efficient and has been used to produce research-grade
results.

At-a-glance:
- given a parametrization of a curve, chunkIE will return
a "chunker" object which stores the description of the
curve into a panel-based format (the curve is discretized into
chunks such that on each chunk a Legendre interpolant in parameter
space is accurate to some prescribed accuracy).
- chunkIE has routines for setting up system matrices
corresponding to logarithmically singular and principal value
type integral equation kernels defined on a chunker
- chunkIE is designed to inter-operate with Ken Ho's fast
linear algebra in MATLAB package (FLAM)
- chunkIE includes various routines for evaluating layer
potentials and functions defined on chunkers
- chunkIE applies a version of Johan Helsing's recursively
compressed inverse preconditioning (RCIP) scheme for effectively
treating problems with corners and multiple junctions.

Detailed documentation is being built [here](https://chunkie.readthedocs.io/en/latest/)

## Installing chunkIE

For other install options, see [the documentation.](https://chunkie.readthedocs.io/en/latest/getchunkie.html)

To compile from git, clone the repository with the submodules 

    git clone --recurse-submodules https://github.com/fastalgorithms/chunkie.git

and run startup.m in the install directory. 
This will download the FLAM and fmm2d submodules, include FLAM in 
the matlab path, and generate the fmm2d mex file if a fortran compiler
exists.

For troubleshooting suggestions, see [the documentation](https://chunkie.readthedocs.io/en/latest/getchunkie.html)

## Using chunkIE

Check out the chunkie/demo folder or the [guide (under construction)](https://chunkie.readthedocs.io/en/latest/guide.html)

## License

chunkIE is copyright 2024 the chunkIE team

chunkIE proper (the contents of the chunkie
folder) is available under the terms of the
BSD 3-clause license, which should have been included
in the distribution (see chunkie/LICENSE.md)

## chunkIE team

chunkIE has benefitted from the contributions of several developers: Travis Askham, 
Manas Rachh, Michael O'Neil, Jeremy Hoskins, Dan Fortunato, Shidong Jiang, 
Fredrik Fryklund, Hai Yang Wang, Hai Zhu, and Tristan Goodwill.

James Bremer and Zydrunas Gimbutas provided generalized Gaussian quadrature rules (chunkie/+chnk/+quadggq)

Many routines were modelled after parts of the legeexps.f library (Copyright Vladimir Rokhlin, Free BSD 3-clause),
FMMLIB2D (Copyright Leslie Greengard and Zydrunas Gimbutas, Free BSD 3-clause), and Johan Helsing's
[RCIP tutorial](https://arxiv.org/abs/1207.6737)

## Citing this software

If you found this software useful, we ask that you please cite the following
works

- {This software} see CITATIONS.cff for details
- {The fast multipole method library} https://github.com/flatironinstitute/fmm2d
- {The fast direct solver library} Ho, Kenneth L. "FLAM: Fast linear algebra in MATLAB-Algorithms for hierarchical matrices." Journal of Open Source Software 5.51 (2020): 1906.
- {Quadrature generation routines} Bremer, James, Zydrunas Gimbutas, and Vladimir Rokhlin. "A nonlinear optimization procedure for generalized Gaussian quadratures." SIAM Journal on Scientific Computing 32.4 (2010): 1761-1788.
- {Corner and multiple junction handling} Helsing, Johan. "Solving integral equations on piecewise smooth boundaries using the RCIP method: a tutorial." Abstract and applied analysis. Vol. 2013. Hindawi, 2013.
- {Fast multipole method paper} Greengard, Leslie, and Vladimir Rokhlin. "A fast algorithm for particle simulations." Journal of computational physics 73.2 (1987): 325-348.

## Contributing

Contributions are welcome. See the issues tab or create
a new issue if there is something you're interested in
bringing to chunkIE. See the
[wiki](https://github.com/fastalgorithms/chunkie/wiki)
for more on the developer process.
