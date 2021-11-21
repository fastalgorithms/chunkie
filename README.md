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

Add the chunkie subfolder to your matlab path.
For some features, you will need to also add the
FLAM library to your path. This is included as
a submodule. You can either clone with the submodules

    git clone --recurse-submodules https://github.com/fastalgorithms/chunkie.git

or initialize the submodules after a git pull.
Once the submodules are initialized, you can
run the setup script (setup.m) from the
chunkie subfolder.

Alternatively, it should work if you already have
a reasonably up-to-date copy of FLAM on your path
see [the FLAM GitHub page](https://github.com/klho/FLAM).
Be sure to add recursively (see MATLAB's genpath
function) to include FLAM's subfolders.


## Using chunkie

Check out the chunkie/demo folder.

## License

chunkie is copyright 2019 the chunkie team

chunkie proper (the contents of the chunkie
folder) is available under the terms of the
BSD 3-clause license, which should have been included
in the distribution (see chunkie/LICENSE.md)

## chunkie team

chunkers:
- Travis Askham
- Manas Rachh
- Michael O'Neil

James Bremer provided generalized Gaussian quadrature rules (chunkie/+chnk/+quad/+brem)

Many routines were modelled after parts of the legeexps.f library (Copyright Vladimir Rokhlin, Free BSD 3-clause) and FMMLIB2D (Copyright Leslie Greengard and Zydrunas Gimbutas, Free BSD 3-clause).

## Contributing

Contributions are welcome. See the issues tab or create
a new issue if there is something you're interested in
bringing to chunkie. See the
[wiki](https://github.com/fastalgorithms/chunkie/wiki)
for more on the developer process.
