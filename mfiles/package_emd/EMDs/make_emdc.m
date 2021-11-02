%MAKE_EMDC  Compiles the C codes for Empirical Mode Decomposition
%
% Note: The compilation can fail on some systems (e.g. MacOS) if Matlab cannot find the C compiler.
% In this case, you should either install a C compiler or check Matlab configuration.
% The configuration files for compilation are mexopts.sh (Unix / MAC OS) and mexopts.bat (Windows)
% use "mex -setup" to choose a configuration file (Unix) or select a compiler (Windows).

function varargout=make_emdc

oldpwd = pwd;
path = fileparts(which('make_emdc'));
cd(path)

if ispc
  cd('src')
end

filelist = {'emdc.c','emdc_fix.c','cemdc.c','cemdc_fix.c','cemdc2.c','cemdc2_fix.c'};

for k = 1:length(filelist)
  file = filelist{k};
  if ispc
    args = {file,'-output', ['../',file(1:end-1),mexext]};
  else
    args = {['src/',file]};
  end
  try
    mex('-DC99_OK',args{:})
    status(k) = 0;
  catch
    try
      mex (args{:})
      status(k) = 1;
    catch
      status(k) = 2;
    end
  end
end

if any(status == 1)
  warning('<complex.h> compiler extension not found. using ANSI C implementation (slower) instead for the following files:')
  for k=find(status==1)
    disp(filelist{k})
  end
end

if any(status == 2)
  warning('The compilation of the following files failed: ')
  for k=find(status==2)
    disp(filelist{k})
  end
end


cd(oldpwd)

if nargout > 0
  varargout = {status};
end