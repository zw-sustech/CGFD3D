function draw_colorbar_fun(pnm_fig,fnmpart,cid,varargin)

flagtick=0; flagtick2=0;

args=varargin; nargs=nargin-3;

if nargs>=1, tickloc=varargin{1};flagtick=1; end
if nargs>=2, tickhalf=varargin{2};flagtick2=1; end

set(gca,'visible','off')

% ----------- south ------------------
dnm='South'; set(cid,'location',dnm); pos=get(cid,'position');

set(cid,'XAxisLocation','bottom');
if flagtick, set(cid,'xtick',tickloc); end
%pause
print(gcf,'-dtiff',[pnm_fig '/' fnmpart 'colorbar_' dnm 'B.tiff']);
set(cid,'XAxisLocation','top');
print(gcf,'-dtiff',[pnm_fig '/' fnmpart 'colorbar_' dnm 'T.tiff']);
set(cid,'xtick',[]);
print(gcf,'-dtiff',[pnm_fig '/' fnmpart 'colorbar_' dnm 'Empty.tiff']);

set(cid,'position',[pos(1),pos(2),pos(3)/2,pos(4)/2]);
if flagtick, set(cid,'xtick',tickloc); end
set(cid,'XAxisLocation','bottom');
print(gcf,'-dtiff',[pnm_fig '/' fnmpart 'colorbar_' dnm 'B2.tiff']);
set(cid,'XAxisLocation','top');
print(gcf,'-dtiff',[pnm_fig '/' fnmpart 'colorbar_' dnm 'T2.tiff']);
set(cid,'xtick',[]);
print(gcf,'-dtiff',[pnm_fig '/' fnmpart 'colorbar_' dnm 'Empty2.tiff']);

set(cid,'position',[pos(1),pos(2),pos(3)/4,pos(4)/4]);
if flagtick, set(cid,'xtick',tickloc); end
if flagtick2, set(cid,'xtick',tickhalf); end

set(cid,'XAxisLocation','bottom');
print(gcf,'-dtiff',[pnm_fig '/' fnmpart 'colorbar_' dnm 'B4.tiff']);
set(cid,'XAxisLocation','top');
print(gcf,'-dtiff',[pnm_fig '/' fnmpart 'colorbar_' dnm 'T4.tiff']);
set(cid,'xtick',[]);
print(gcf,'-dtiff',[pnm_fig '/' fnmpart 'colorbar_' dnm 'Empty4.tiff']);

if flagtick
   set(cid,'xtick',[tickloc(1),tickloc(end)]);
else
   tik=get(cid,'xlim');
   if tik(1)<0 & tik(end)>0
      tik=[tik(1),0,tik(end)];
   else
      tik=[tik(1),tik(end)]; 
   end
   set(cid,'xtick',tik);
end
set(cid,'XAxisLocation','bottom');
print(gcf,'-dtiff',[pnm_fig '/' fnmpart 'colorbar_' dnm 'B43.tiff']);
set(cid,'XAxisLocation','top');
print(gcf,'-dtiff',[pnm_fig '/' fnmpart 'colorbar_' dnm 'T43.tiff']);

% ----------- east ------------------
dnm='East'; set(cid,'location',dnm); pos=get(cid,'position');

set(cid,'YAxisLocation','left');
if flagtick, set(cid,'ytick',tickloc); end
print(gcf,'-dtiff',[pnm_fig '/' fnmpart 'colorbar_' dnm 'L.tiff']);
set(cid,'YAxisLocation','right');
print(gcf,'-dtiff',[pnm_fig '/' fnmpart 'colorbar_' dnm 'R.tiff']);
ytk=get(cid,'ytick'); set(cid,'ytick',[]);
print(gcf,'-dtiff',[pnm_fig '/' fnmpart 'colorbar_' dnm 'Empty.tiff']);
set(cid,'ytick',ytk);

set(cid,'position',[pos(1),pos(2),pos(3)/2,pos(4)/2]);
if flagtick, set(cid,'ytick',tickloc); end
set(cid,'YAxisLocation','left');
print(gcf,'-dtiff',[pnm_fig '/' fnmpart 'colorbar_' dnm 'L2.tiff']);
set(cid,'YAxisLocation','right');
print(gcf,'-dtiff',[pnm_fig '/' fnmpart 'colorbar_' dnm 'R2.tiff']);
ytk=get(cid,'ytick'); set(cid,'ytick',[]);
print(gcf,'-dtiff',[pnm_fig '/' fnmpart 'colorbar_' dnm 'Empty2.tiff']);
set(cid,'ytick',ytk);

set(cid,'position',[pos(1),pos(2),pos(3)/4,pos(4)/4]);
if flagtick, set(cid,'ytick',tickloc); end
if flagtick2, set(cid,'ytick',tickhalf); end

set(cid,'YAxisLocation','left');
print(gcf,'-dtiff',[pnm_fig '/' fnmpart 'colorbar_' dnm 'L4.tiff']);
set(cid,'YAxisLocation','right');
print(gcf,'-dtiff',[pnm_fig '/' fnmpart 'colorbar_' dnm 'R4.tiff']);
ytk=get(cid,'ytick'); set(cid,'ytick',[]);
print(gcf,'-dtiff',[pnm_fig '/' fnmpart 'colorbar_' dnm 'Empty4.tiff']);
set(cid,'ytick',ytk);

if flagtick
   set(cid,'ytick',[tickloc(1),tickloc(end)]);
else
   tik=get(cid,'ylim');
   if tik(1)<0 & tik(end)>0
      tik=[tik(1),0,tik(end)];
   else
      tik=[tik(1),tik(end)]; 
   end
   set(cid,'ytick',tik);
end
set(cid,'YAxisLocation','left');
print(gcf,'-dtiff',[pnm_fig '/' fnmpart 'colorbar_' dnm 'L43.tiff']);
set(cid,'YAxisLocation','right');
print(gcf,'-dtiff',[pnm_fig '/' fnmpart 'colorbar_' dnm 'R43.tiff']);

