# chunkie: a MATLAB integral equation toolbox

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
linear algebra in MATLAB package (FLAM)
- chunkie includes various routines for evaluating layer
potentials and functions defined on chunkers

## Installing chunkie

Clone the repository with the submodules 

    git clone --recurse-submodules https://github.com/fastalgorithms/chunkie.git

and run startup.m in the install directory. 
This will download the FLAM and fmm2d submodules, include FLAM in 
the matlab path, and generate the fmm2d mex file if a fortran compiler
exists. 

## Using chunkie

Check out the chunkie/demo folder.

## License

chunkie is copyright 2019 the chunkie team

chunkie proper (the contents of the chunkie
folder) is available under the terms of the
BSD 3-clause license, which should have been included
in the distribution (see chunkie/LICENSE.md)

## Installation notes

- The fmm2d mex installation is currently not supported on Windows, to
  complete the mex installation, follow instructions on the [fmm2d documentation](https://fmm2d.readthedocs.io/en/latest/install.html) 
- fmm2d mex installation depends on gfortran. In case a compiler is not
  found, the installation will be skipped. To obtain gfortran on MacOS,
  get Xcode, and command line tools using

    xcode-select --install
  Then install Homebrew from https://brew.sh, and finally install
  gfortran using
    
    brew install gcc

  On ubuntu linux run

    sudo apt-get install make build-essential gfortran

  On fedora/centOS linux run

    sudo yum install make gcc gcc-c++ gcc-gfortran libgomp

- If installing without submodules, chunkie depends on [FLAM](https://github.com/klho/FLAM), 
and optionally on the
[fmm2d](https://github.com/flatironinstitute/fmm2d) repository. Parts of
the library will not function without FLAM and its subdirectories included in the matlab path.


## chunkie team

chunkers:
- Travis Askham
- Manas Rachh
- Michael O'Neil
- Jeremy Hoskins
- Dan Fortunato

James Bremer provided generalized Gaussian quadrature rules (chunkie/+chnk/+quad/+brem)

Many routines were modelled after parts of the legeexps.f library (Copyright Vladimir Rokhlin, Free BSD 3-clause) and FMMLIB2D (Copyright Leslie Greengard and Zydrunas Gimbutas, Free BSD 3-clause).

## Contributing

Contributions are welcome. See the issues tab or create
a new issue if there is something you're interested in
bringing to chunkie. See the
[wiki](https://github.com/fastalgorithms/chunkie/wiki)
for more on the developer process.
