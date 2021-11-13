%INSTALL_EMD.M install the EMD package
% add the following locations to Matlab's path:
% <EMD_BASEDIR>
% <EMD_BASEDIR>/emd
% <EMD_BASEDIR>/utils
% <EMD_BASEDIR>/examples/NSIP2003
% <EMD_BASEDIR>/examples/SPL2007
%
% where <EMD_BASEDIR> is the directory containing install_emd.m
%
% and compiles the C codes
%
% IMPORTANT: After running INSTALL_EMD you must run the "savepath" command to save the installation
% but be careful that if you previously removed parts of the path (using e.g. the "rmpath" command) 
% these will be permanently removed after you run "savepath"

function install_emd

base_dir = fileparts(which('install_emd'));
addpath(base_dir)
addpath([base_dir,'/EMDs'])
addpath([base_dir,'/utils'])
addpath([base_dir,'/examples/NSIP2003'])
addpath([base_dir,'/examples/SPL2007'])

status=make_emdc;
if all(status==0)
  fprintf('\nCompilation successfull\n\n')
elseif all(status~=2)
  fprintf('\nCompilation successfull.\n\n')
  fprintf('Some codes can run faster if they are compiled with a C compiler\n')
  fprintf('that handles the C99 complex data type ("complex.h"). See details above.\n\n')
else
  fprintf('\nSome errors occurred during compilation. See details above.\n\n')
end

disp('Installation complete. Run index_emd for a list of functions.')
fprintf('\n')
disp('IMPORTANT: After running INSTALL_EMD you must run the "savepath" command to save the installation')
disp('but be careful that if you previously removed parts of the path (using e.g. the "rmpath" command)')
disp('these will be permanently removed after you run "savepath"')
