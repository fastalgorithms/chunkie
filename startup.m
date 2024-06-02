function startup(opts)
%STARTUP chunkIE startup routine 
%
% adds necessary folders to path. checks if FLAM and fmm2d directories are
% in the right locations. installs fmm2d if needed and possible
%
% Optional input:
%
% opts = options structure (defaults)
%   opts.testfmm - boolean (false), test fmm intallation
%   opts.fmmremex - boolean (false), force recompile of mex routine
%   opts.fmmrecompile - boolean (false), force recompile of fmm2d library
%

if nargin < 1
    opts = [];
end

fmmremex = false;
if isfield(opts,'fmmremex')
    fmmremex = opts.fmmremex;
end

fmmrecompile = false;
if isfield(opts,'fmmrecompile')
    fmmrecompile = opts.fmmrecompile;
end

fmmremex = or(fmmremex,fmmrecompile);

testfmm = false;
if isfield(opts, 'testfmm')
  testfmm = opts.testfmm;
end

% add paths 
addpath('./chunkie')
if(exist('./chunkie/FLAM/startup.m','file'))
    run 'chunkie/FLAM/startup.m';
else
    msg = "STARTUP: warning FLAM not found in usual location " + ...
        "Check that the submodule was included. ";
    warning(msg);
end

% install fmm2d if needed, add to path if found
if(exist('chunkie/fmm2d/matlab','dir'))
    cd './chunkie/fmm2d';
    addpath './matlab';
    icheck = exist(['fmm2d.' mexext], 'file');
    if icheck ~=3 || fmmremex || fmmrecompile
        if ismac || isunix
            [status,cmdout] = system('which gfortran');
            if(~status)
                fprintf('fortran compiler found at: %s\n',cmdout);
                iffmm = true;
                path1 = getenv('PATH');
                cmdout2 = extractBefore(cmdout, 'gfortran');
                path1 = [path1 cmdout2];
                setenv('PATH', path1);
                if ismac
                    !cp -f chunkie/fmm2d/make.inc.macos.gnu.x86-64 chunkie/fmm2d/make.inc;
                end
                if fmmrecompile
                    !make clean
                end
                !make matlab;  
                testfmm = true; % test if new install

            else    
                msg = "STARTUP: unable to find a suitable compiler for FMM2D. " + ...
                    "See manual install instructions on github";
                warning(msg)
            end
                
        else
            msg = "STARTUP: automatic install not supported on your operating system. " + ...
                "See manual install instructions on github";
            warning(msg)
        end        
    end
    cd matlab;
    if testfmm
        runtests;
    end
    cd ../../../;
else
    if fmmremex || fmmrecompile
        msg = "STARTUP: fmm2d reinstall requested but source files not found in usual location" + ...
        "Check that the submodule was included. ";
        warning(msg);
    end
    dir1 = pwd();
    icheck = exist(['fmm2d.' mexext], 'file');
    if icheck == 3
        fname = which(['fmm2d.' mexext]);
        dir0 = dir(fname);
        msg = "STARTUP: fmm2d mex file found on PATH but not in usual location " + ...
        "adding " + dir0.folder + " to path";
        warning(msg);
        addpath(dir0.folder);
        
        if testfmm
            msg = "STARTUP: testing fmm2d not supported if not installed in usual location";
            warning(msg);
        end
    else
        msg = "STARTUP: warning fmm2d folder not found in usual location " + ...
            "and no fmm2d mex file is known to MATLAB. If mex file exists " + ...
            " in another location, add that location to path. " + ...
        "Else, check that the submodule was included. " + ...
        "Else, check that startup is run from its folder. ";
        
        warning(msg);
    end
end
