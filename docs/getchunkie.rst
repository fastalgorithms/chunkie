.. role:: matlab(code)
   :language: matlab   

Get chunkIE
============

Installation From Source 
---------------------------

Installation directly from source is not recommended for Windows
machines. See source code with binaries below.

- Get the source code:
  
  * You can download the latest `stable release <https://github.com/fastalgorithms/chunkie/releases/download/v1.0.0/chunkie-v1.0.0.zip>`_
    
  * For the latest and greatest, you can download from the git repository using
  
  .. code:: bash
	  
     git clone --recurse-submodules https://github.com/fastalgorithms/chunkie.git

- Once you have the source code, you run startup.m in the install directory. This will compile fmm2d binaries and create a mex file, which are the main source of installation trouble (see troubleshooting below). Installation succeeded if the fmm tests run successfully (startup will automatically run these after installing).
- The chunkie library will function without the fmm2d binaries, but it will
  be slower for some calculations. 
  
Install with Precompiled Windows Binaries
------------------------------------------

Windows compilation of fmm2d binaries and mex files has many issues.
For intel systems, we have zip files with `pre-compiled binaries Windows
available <https://github.com/fastalgorithms/chunkie/releases/tag/v1.0.0>`_.

Instructions:

- Select the zip file whose name has MATLAB version 2020 if your version is 2022 or older. Select 2023 otherwise.

- The x86 version should work on most intel machines. The avx2 specification is also commonly available and will be faster if it is compatible with your machine.

- Download the file and unzip. We recommend testing the fmm binaries by running the following in the top-level chunkie directory:

.. code:: matlab

   startup(struct("testfmm",true))


Troubleshooting
-----------------

- "{chunkIE function name} is not found".
  MATLAB needs to know the location of the chunkIE source files to be able
  to use chunkIE routines. You can add all of the necessary files to path by
  running startup in the top-level chunkie directory. There are ways to
  `add these locations to the MATLAB path on startup <https://www.mathworks.com/help/matlab/matlab_env/add-folders-to-matlab-search-path-at-startup.html>`_
- Windows install issues. The fmm2d mex installation is currently not supported on Windows. Either
  use one of the pre-compiled options above or attempt to
  complete the mex installation by follow instructions on the `fmm2d documentation <https://fmm2d.readthedocs.io/en/latest/install.html>`_
- "CHUNKIE STARTUP: unable to find a suitable compiler for FMM2D."
  fmm2d mex installation depends on gfortran. In case a compiler is not
  found, the installation will be skipped. To install dependencies follow the procedure below based on your OS
- "mex" redirects to LaTeX compilation. This can happen on either MacOS or Linux, and in both situations the path
  to mex is not correctly set. We recommend doing a manual install in the fmm2d folder. 
  Find your MATLAB installation directory using ``${MATLABROOT}`` on linux, 
  and on MacOS, typically MATLAB is installed at ``/Applications/MATLAB_R(version number)/``.
  The mex executable can then be found in the ``bin`` directory. Navigate to chunkie/fmm2d. 
  Create a make.inc file by copying one of the template files (make.inc.*) into make.inc. For example,
  on MacOS you might do 

  .. code:: bash

        cp make.inc.macos.gnu make.inc

  There are other make.inc templates for special situations; see the fmm2d folder. 
  Then add the following to a new line at the end of make.inc: ``MEX=<path to mex binary>``. Note: you 
  should use the absolute path, not a relative path. Then, run 

  .. code:: bash
  
        make matlab


- "CHUNKIE STARTUP: unable to find a suitable compiler for FMM2D". You may not have gfortran installed (if you do and you're on MacOS, see further below). To install gfortran, try:

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

- "CHUNKIE STARTUP: unable to find a suitable compiler for FMM2D". If you have gfortran installed but you still get this message on MacOS, try starting MATLAB from the command line (as opposed to opening MATLAB via spotlight). To do this, you may have to find the location of the MATLAB binary in the bin folder of MATLAB's install directory and run that path directly. See also `instructions from MathWorks <https://www.mathworks.com/help/matlab/ref/matlabmacos.html>`_

  .. code:: bash
    
        ./<path to MATLAB executable>

- "{hypoct, rskelf, ifmm} is not found". These are part of the FLAM library. If you've already run startup.m, then it may be that you downloaded from git but forgot to recurse submodules. Do the download from git again and be sure to include the submodules.
