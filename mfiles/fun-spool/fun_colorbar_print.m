function cid=fun_colorbar_print(sid,FIG_ROOT,fnmpart,dfmt,filefmt,ticks,clim)

if ~isdir(FIG_ROOT), mkdir(FIG_ROOT); end

drvfmt=['-d' dfmt];
if exist('filefmt','var')
   fnmfmt=filefmt;
else
   fnmfmt=dfmt;
end

colorbar off
cid=colorbar('vert','location','SouthOutSide');
set(gca,'visible','off');
set(sid,'visible','off');

if exist('ticks','var')
   set(cid,'xtick',ticks);
end
scl_tick=get(cid,'xtick');

% -- print sourth --
if exist('clim','var')
   set(cid,'xlim',clim);
end
scl_clim=get(cid,'xlim'); v0=scl_clim(2);
scl_tick(find(scl_tick>v0))=[];
if v0-scl_tick(end)<(scl_tick(end)-scl_tick(end-1))/2, scl_tick(end)=[]; end
scl_tick=[scl_tick v0];
scl_tick2=fliplr(scl_tick(end:-2:1));
scl_tick4=fliplr(scl_tick(end:-4:1));

print_south([]);
set(gcf,'InvertHardcopy','off')
print_south('wbg_');
set(gcf,'InvertHardcopy','on')
% -- restore display --
set(gca,'visible','on');
set(sid,'visible','on');

% -- print east --
colorbar off
cid=colorbar('vert','location','EastOutSide');
set(gca,'visible','off');
set(sid,'visible','off');

if exist('clim','var')
   set(cid,'ylim',clim);
end
print_east([]);
set(gcf,'InvertHardcopy','off')
print_east('wbg_');
set(gcf,'InvertHardcopy','on')

% -- restore display --
set(gca,'visible','on');
set(sid,'visible','on');

% ----------- south ------------------
function print_south(wbg)

dnm='South'; set(cid,'location',dnm); pos=get(cid,'position');
   
% whole
if exist('ticks','var'), set(cid,'xtick',scl_tick); end

set(cid,'XAxisLocation','bottom');
drawnow
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_colorbar_' wbg dnm 'B.' fnmfmt]);
drawnow
set(cid,'XAxisLocation','top');
drawnow
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_colorbar_' wbg dnm 'T.' fnmfmt]);
drawnow
   
% 1/2
set(cid,'position',[pos(1),pos(2),pos(3)/2,pos(4)/2]);
if exist('ticks','var'), set(cid,'xtick',scl_tick2); end

set(cid,'XAxisLocation','bottom');
drawnow
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_colorbar_' wbg dnm 'B2.' fnmfmt]);
drawnow
set(cid,'XAxisLocation','top');
drawnow
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_colorbar_' wbg dnm 'T2.' fnmfmt]);
drawnow
   
% 1/4
set(cid,'position',[pos(1),pos(2),pos(3)/4,pos(4)/4]);
if exist('ticks','var'), set(cid,'xtick',scl_tick4); end

set(cid,'XAxisLocation','bottom');
drawnow
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_colorbar_' wbg dnm 'B4.' fnmfmt]);
drawnow
set(cid,'XAxisLocation','top');
drawnow
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_colorbar_' wbg dnm 'T4.' fnmfmt]);
drawnow
   
% 1/4+0
set(cid,'xtick',scl_clim);
set(cid,'XAxisLocation','bottom');
drawnow
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_colorbar_' wbg dnm 'B43.' fnmfmt]);
drawnow
set(cid,'XAxisLocation','top');
drawnow
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_colorbar_' wbg dnm 'T43.' fnmfmt]);
drawnow
   
% 1/4+notick
set(cid,'xtick',[]);
drawnow
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_colorbar_' wbg dnm 'E4.' fnmfmt]);
drawnow
end

% ----------- east ------------------
function print_east(wbg)

dnm='East'; set(cid,'location',dnm); pos=get(cid,'position');
   
% whole
if exist('ticks','var'), set(cid,'ytick',scl_tick); end

set(cid,'YAxisLocation','left');
drawnow
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_colorbar_' wbg dnm 'L.' fnmfmt]);
drawnow
set(cid,'YAxisLocation','right');
drawnow
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_colorbar_' wbg dnm 'R.' fnmfmt]);
drawnow
   
% 1/2
set(cid,'position',[pos(1),pos(2),pos(3)/2,pos(4)/2]);
if exist('ticks','var'), set(cid,'ytick',scl_tick2); end

set(cid,'YAxisLocation','left');
drawnow
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_colorbar_' wbg dnm 'L2.' fnmfmt]);
drawnow
set(cid,'YAxisLocation','right');
drawnow
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_colorbar_' wbg dnm 'R2.' fnmfmt]);
drawnow
   
% 1/4
set(cid,'position',[pos(1),pos(2),pos(3)/4,pos(4)/4]);
if exist('ticks','var'), set(cid,'ytick',scl_tick4); end

set(cid,'YAxisLocation','left');
drawnow
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_colorbar_' wbg dnm 'L4.' fnmfmt]);
drawnow
set(cid,'YAxisLocation','right');
drawnow
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_colorbar_' wbg dnm 'R4.' fnmfmt]);
drawnow
   
% 1/4+0
set(cid,'ytick',scl_clim);
set(cid,'YAxisLocation','left');
drawnow
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_colorbar_' wbg dnm 'L43.' fnmfmt]);
drawnow
set(cid,'YAxisLocation','right');
drawnow
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_colorbar_' wbg dnm 'R43.' fnmfmt]);
drawnow
   
% 1/4+notick
set(cid,'ytick',[]);
drawnow
print(gcf,drvfmt,[FIG_ROOT '/' fnmpart '_colorbar_' wbg dnm 'E4.' fnmfmt]);
drawnow
end

end
