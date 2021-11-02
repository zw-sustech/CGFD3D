function fun_print(sid,FIG_ROOT,fnmpart,dfmt,filefmt)

if ~isdir(FIG_ROOT), mkdir(FIG_ROOT); end

drvfmt=['-d' dfmt];
if exist('filefmt','var')
   fnmfmt=filefmt;
else
   fnmfmt=dfmt;
end

% ----------------------- print -----------------------

print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '.' fnmfmt]);

set(gca,'visible','off')
set(sid,'visible','on')
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_objonly.' fnmfmt]);

set(gca,'visible','on')
set(sid,'visible','off')
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_gcaonly.' fnmfmt]);

%oldstat=get(gcf,'InvertHardcopy');
%oldgcastat=get(gca,'visible');

set(gcf,'InvertHardcopy','off')

set(gca,'visible','on')
set(sid,'visible','on')
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_wbgc.' fnmfmt]);

set(gca,'visible','off')
set(sid,'visible','on')
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_wbgc_objonly.' fnmfmt]);

set(gca,'visible','on')
set(sid,'visible','off')
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_wbgc_gcaonly.' fnmfmt]);

set(gcf,'InvertHardcopy','on')
