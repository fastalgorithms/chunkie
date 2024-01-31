function [] = startup(varargin)

opts = [];
if(nargin == 1)
   opts = nargin;
end

ifflam = true;
if(isfield(opts,'ifflam'))
   ifflam = opts.ifflam;
end


addpath('./chunkie')
if(ifflam)
  if(exist('chunkie/FLAM/startup.m'))
    run 'chunkie/FLAM/startup.m';
  end
end

iffmm0 = true;


if(isfield(opts,'iffmm'))
  iffmm0 = opts.iffmm;
end

iffmm = true;
if(iffmm0)
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
        else
            fprintf('Fortran compiler not found\n');
            iffmm = false;
        end

    else
       fprintf('Fortran installations not supported through startup.m\n');
       fprintf('Follow manual installation instructions on github\n');
       iffmm = false;
    end
end
        



if(iffmm)
  if(exist('chunkie/fmm2d/matlab','dir'))
      cd './chunkie/fmm2d';
      !make clean;
      !make matlab;  
      addpath './matlab';
      cd matlab;
      runtests;
      cd ../../../;
  else
      fprintf('Fmm installation not in standard location\n');
      fprintf('Manual installation required');
  end
end
