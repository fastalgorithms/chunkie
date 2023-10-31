Get Chunkie
============

Installation Instructions
--------------------------

Clone the repository with the submodules 

.. code:: bash
	  
    git clone --recurse-submodules https://github.com/fastalgorithms/chunkie.git

and run startup.m in the install directory. 
This will download the FLAM and fmm2d submodules, include FLAM in 
the matlab path, and generate the fmm2d mex file if a fortran compiler
exists. 


Troubleshooting
-----------------

- The fmm2d mex installation is currently not supported on Windows, to
  complete the mex installation, follow instructions on the `fmm2d documentation <https://fmm2d.readthedocs.io/en/latest/install.html>`_
- fmm2d mex installation depends on gfortran. In case a compiler is not
  found, the installation will be skipped. To install dependencies follow the procedure below based on your OS
  
  * MacOS
  
    Get xcode, command line tools by running

    .. code:: bash
    
        xcode-select --install
    
    Then install Homebrew from https://brew.sh, and finally install gfortran using

    .. code:: bash
  
        brew install gcc

  * Ubuntu linux

    .. code:: bash

       sudo apt-get install make build-essential gfortran

  * Fedora/centOS linux

    .. code:: bash
    
       sudo yum install make gcc gcc-c++ gcc-gfortran libgomp

- If installing without submodules, chunkie depends on `FLAM <https://github.com/klho/FLAM>`_, and optionally on the `fmm2d <https://github.com/flatironinstitute/fmm2d>`_ repository. Parts of the library will not function without FLAM and its subdirectories included in the matlab path.
