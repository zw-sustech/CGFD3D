% UNINSTALL_EMD.M uninstall the EMD package
% reverts the modifications made by INSTALL_EMD:
% removes the following locations from Matlab's path:
% <EMD_BASEDIR>
% <EMD_BASEDIR>/utils
% <EMD_BASEDIR>/examples/NSIP2003
% <EMD_BASEDIR>/examples/SPL2007
%
% where <EMD_BASEDIR> is the directory containing uninstall_emd.m
%
% UNINSTALL_EMD also interactively proposes to remove the files from the hard drive
%
% IMPORTANT: After running UNINSTALL_EMD you must run the "savepath" command to save the modifications made 
% to Matlab's path but be careful that if you previously removed parts of the path (using e.g. the "rmpath" command) 
% these will be permanently removed after you run "savepath"
function uninstall_emd
base_dir = fileparts(which('uninstall_emd'));
rmpath(base_dir)
rmpath([base_dir,'/EMDs'])
rmpath([base_dir,'/utils'])
rmpath([base_dir,'/examples/NSIP2003'])
rmpath([base_dir,'/examples/SPL2007'])

reply = input('Would you like to remove the files from the hard drive? Y/N [N]: ', 's');
if isempty(reply)
  reply = 'n';
end

if ~any(strcmpi(reply,{'n','y'}))
  disp(['Unsupported answer ',reply,'. Doing nothing.'])
end

if strcmpi(reply,'y')
  dir_list = {};
  file_list = textread([base_dir,'/ls-R'],'%s');
  for file = file_list'
    file = char(file);
    if file(1) == '.'
      sub_dir = file(2:end-1);
      dir_list{end+1} = [base_dir,sub_dir];
    else
      filename = [base_dir,sub_dir,'/',file];
      if ~isdir(filename)
        [p,n,e] = fileparts(filename);
        if strcmp(e,'mexglx')
          filename = [p,n,mexext];
        end
        if exist(filename,'file')
          delete(filename)
        else
          warning(['File not found: ',filename])
        end
      end
    end
  end
  for dir = fliplr(dir_list)
    dir = char(dir);
    [status,message] = rmdir(dir);
  end
end

disp('IMPORTANT: After running UNINSTALL_EMD you must run the "savepath" command to save the modifications made')
disp('to Matlab''s path but be careful that if you previously removed parts of the path (using e.g. the "rmpath" command)')
disp('these will be permanently removed after you run "savepath"')