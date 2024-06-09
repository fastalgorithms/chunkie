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
    msg = "CHUNKIE STARTUP: warning FLAM not found in usual location\n " + ...
        "Check that the submodule was included. ";
    warning('CHUNKIESTARTUP:flamnotfoundwarn',msg);
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
                fprintf('------- chunkIE startup: building fmm2d ----- \n');
                fprintf('fortran compiler found at: %s\n',cmdout);
                iffmm = true;
                path1 = getenv('PATH');
                cmdout2 = extractBefore(cmdout, 'gfortran');
                path1 = [path1 cmdout2];
                setenv('PATH', path1);
                if ismac
                    [~, result] = system('uname -m');
                    if strcmpi(result, 'x86_64')
                        !cp -f make.inc.macos.gnu make.inc;
                    else
                        !cp -f make.inc.macos_arm64.gnu make.inc;
                    end
                end
                if fmmrecompile
                    !make clean
                end
                !make matlab;  
                testfmm = true; % test if new install

            else    
                msg = "CHUNKIE STARTUP: unable to find a suitable compiler for FMM2D. " + ...
                    "See manual install instructions on github";
                warning('CHUNKIESTARTUP:compilerwarn',msg)
            end
                
        else
            msg = "CHUNKIE STARTUP: automatic install not supported on your operating system. " + ...
                "See manual install instructions on github";
            warning('CHUNKIESTARTUP:OSnosupport',msg)
        end        
    end
    cd matlab;
    if testfmm
        runtests;
    end
    cd ../../../;
else
    if fmmremex || fmmrecompile
        msg = "CHUNKIE STARTUP: fmm2d reinstall requested but source files not found in usual location.\n" + ...
        "Check that the submodule was included. ";
        warning('CHUNKIESTARTUP:fmm2dsourcenotfoundwarn',msg);
    end
    icheck = exist(['fmm2d.' mexext], 'file');
    if icheck == 3
        msg = "CHUNKIE STARTUP: fmm2d mex file found on PATH but not in usual location " + ...
        "There may be unexpected behavior if this is a different fmm2d version or if other fmm2d MATLAB files are not known to MATLAB.";
        warning('CHUNKIESTARTUP:fmm2ddifflocwarn',msg);
        
        if testfmm
            msg = "CHUNKIE STARTUP: testing fmm2d not supported if not installed in usual location";
            warning('CHUNKIESTARTUP:fmm2dnotestwarn',msg);
        end
    else
        msg = ['CHUNKIE STARTUP: warning fmm2d folder not found in usual location ' ...
            'and no fmm2d mex file is known to MATLAB.\n If mex file exists ' ...
            'in another location, add that location to path.\n '...
            'Else, check that the submodule was included. '];
        
        warning('CHUNKIESTARTUP:fmm2dmexnotfoundwarn',msg);
    end
end
