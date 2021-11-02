%P2   plotting utility for SAC files
%
%     plot SAC seismograms read into matlab with rsac.m on top
%     of each other
% 
%     usage:  p2(file1, file2, ... , fileN, [time1 time2])
%
%     where file1,file2,... represents any matlab variable
%     that is a seismogram read in with rsac.m
%
%     [time1 time2] are optional.  If these are not specified
%     the entire time trace will be displayed, otherwise
%     all traces will be plotted over this time window.
%
%     The utility draws a legend where the traces are labeled by
%     the matlab variable name.  The legend box can be moved by
%     left-clicking on the box, and dragging it to the desired
%     location.  The legend labels can also be edited by double-
%     clicking on them.
%
%     Example:
%
%     To plot aak and casy0 for the entire time series on top of
%     each other:
%
%     p2(aak,casy0) 
%
%     To plot aak, casy0, and hrv for the time range of 0 to 300
%     seconds on top of each other:
%
%     p2(aak,casy0,hrv,[0 300]) 
%
%     by Michael Thorne (5/2004)   mthorne@asu.edu

function p2(varargin) 

set(gcf,'Name','P2 -- SAC Seismogram Plotting Utility', ...
    'NumberTitle','off','Color',[.8 .8 .8], ...
    'Pointer','crosshair','PaperOrientation','landscape', ...
    'PaperPosition',[.5 2 10 4.5],'PaperType','usletter');

[a,b]=size(varargin{nargin});
junk=0;
loopend=nargin;
if a==1 & b==2
  limits=varargin{nargin};
  xaxmin=limits(1,1);
  xaxmax=limits(1,2); 
  junk=1;
  loopend=nargin-1;
end

colormap(lines)
pcolors=colormap;

ymin=9999999;
ymax=-9999999;

for i=1:loopend
  file=varargin{i};

    c1=pcolors(i,1); c2=pcolors(i,2); c3=pcolors(i,3);
    plot(file(:,1),file(:,2),'Color',[c1 c2 c3])
    set(gca,'Xcolor',[.1 .1 .1],'Ycolor',[.1 .1 .1], ...
        'FontName','times','FontWeight','light', ...
        'TickDir','out','FontSize',8)
    grid on

    if junk == 1
      bb=find(file(:,1)<xaxmax & file(:,1)>xaxmin);
      temp=file(bb,2);
      minny=min(temp);
      maxxy=max(temp);
      if minny < ymin
        ymin = minny;
      end
      if maxxy > ymax
        ymax = maxxy;
      end
    end

    inames{i}=inputname(i);

    hold on
end

if junk == 1
  axis([xaxmin xaxmax ymin ymax])
end

xlabel('Time (sec)')
legend(inames,0)

hold off
