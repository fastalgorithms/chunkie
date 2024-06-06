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
corresponding to logarithmically singular integral equation
kernels defined on a chunker
- chunkIE is designed to inter-operate with Ken Ho's fast
linear algebra in MATLAB package (FLAM)
- chunkIE includes various routines for evaluating layer
potentials and functions defined on chunkers

## Installing chunkIE

Clone the repository with the submodules 

    git clone --recurse-submodules https://github.com/fastalgorithms/chunkie.git

and run startup.m in the install directory. 
This will download the FLAM and fmm2d submodules, include FLAM in 
the matlab path, and generate the fmm2d mex file if a fortran compiler
exists. 

## Using chunkIE

Check out the chunkie/demo folder.

## License

chunkIE is copyright 2019 the chunkIE team

chunkIE proper (the contents of the chunkie
folder) is available under the terms of the
BSD 3-clause license, which should have been included
in the distribution (see chunkie/LICENSE.md)

## Installation notes

- The fmm2d mex installation is currently not supported on Windows, to
  complete the mex installation, follow instructions on the [fmm2d documentation](https://fmm2d.readthedocs.io/en/latest/install.html) 
- fmm2d mex installation depends on gfortran. In case a compiler is not
  found, the installation will be skipped. To install dependencies follow the procedure below based on your OS
  
  * MacOS
  
    Get xcode, command line tools by running
    
        xcode-select --install
    
    Then install Homebrew from https://brew.sh, and finally install gfortran using
  
        brew install gcc

  * Ubuntu linux

        sudo apt-get install make build-essential gfortran

  * Fedora/centOS linux

        sudo yum install make gcc gcc-c++ gcc-gfortran libgomp

- If installing without submodules, chunkIE depends on [FLAM](https://github.com/klho/FLAM), 
and optionally on the
[fmm2d](https://github.com/flatironinstitute/fmm2d) repository. Parts of
the library will not function without FLAM and its subdirectories included in the matlab path.


## chunkIE team

chunkIE has benefitted from the contributions of several developers: Travis Askham, 
Manas Rachh, Michael O'Neil, Jeremy Hoskins, Dan Fortunato, Shidong Jiang, 
Fredrik Fryklund, Hai Yang Wang, Hai Zhu, and Tristan Goodwill.

James Bremer and Zydrunas Gimbutas provided generalized Gaussian quadrature rules (chunkie/+chnk/+quadggq)

Many routines were modelled after parts of the legeexps.f library (Copyright Vladimir Rokhlin, Free BSD 3-clause),
FMMLIB2D (Copyright Leslie Greengard and Zydrunas Gimbutas, Free BSD 3-clause), and Johan Helsing's
RCIP tutorial (https://arxiv.org/abs/1207.6737)

## Citing this software

If you found this software useful, we ask that you please cite the following
works

- [This software] see CITATIONS.cff for details
- [The fast multipole method library] https://github.com/flatironinstitute/fmm2d
- [The fast direct solver library] Ho, Kenneth L. "FLAM: Fast linear algebra in MATLAB-Algorithms for hierarchical matrices." Journal of Open Source Software 5.51 (2020): 1906.
- [Quadrature generation routines] Bremer, James, Zydrunas Gimbutas, and Vladimir Rokhlin. "A nonlinear optimization procedure for generalized Gaussian quadratures." SIAM Journal on Scientific Computing 32.4 (2010): 1761-1788.
- [Corner and multiple junction handling] Helsing, Johan. "Solving integral equations on piecewise smooth boundaries using the RCIP method: a tutorial." Abstract and applied analysis. Vol. 2013. Hindawi, 2013.
- [Fast multipole method paper] Greengard, Leslie, and Vladimir Rokhlin. "A fast algorithm for particle simulations." Journal of computational physics 73.2 (1987): 325-348.

## Contributing

Contributions are welcome. See the issues tab or create
a new issue if there is something you're interested in
bringing to chunkIE. See the
[wiki](https://github.com/fastalgorithms/chunkie/wiki)
for more on the developer process.
