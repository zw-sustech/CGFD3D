function fid = func_figure( wid, hwratio, render )

%=====================================================
% figure setting
%=====================================================
if wid == 1
    Winch=6.83; %- two column figure (AGU)
elseif wid == 0.5
    Winch=3.33;
elseif wid == 0.75
    Winch=3.33*1.5;
else
    Winch=6.83; %- default two column figure (AGU)
end

Hinch = Winch * hwratio;
if Hinch > 9.66
    disp(['Height=',num2str(Hinch),' exceeds the limit, reset to 9.66']);
    Hinch=9.66; %- largest dep (AGU)
end

Wline=0.5;

fid = figure;
%=== painters for eps file, zbuffer for bitmap file ===
if exist('render','var')
    set(gcf,'renderer',render);
else
    set(gcf,'renderer','painters');
end

%=== paper size ===
%set(gcf,'PaperPositionMode','auto') % the same size as on screen and centered
set(gcf,'PaperPositionMode','manual') % honors the PaperPosition
set(gcf, 'PaperUnits', 'inches')
set(gcf,'PaperSize',[Winch Hinch])    
%=== printed figure to lower left corner of the page and width, height ===
set(gcf,'PaperPosition',[0 0 Winch Hinch])
%=== figure size on screen ===
set(gcf,'Units','inches')
set(gcf,'Position',[2 2 Winch Hinch])
%=== retain the background color ===
%set(gcf,'color','blue')  % bg color of figure
%set(gca,'color','red')   % bg color of axes
%set(gcf,'InvertHardCopy','off')

set(gcf,'defaultaxesfontname','Times New Roman')
set(gcf,'defaultaxesfontsize',8);
set(gcf,'defaulttextfontsize',8);
set(gcf,'defaultaxeslinewidth',0.5);
set(gcf,'defaultlinelinewidth',0.5);


