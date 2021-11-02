%P1   plotting utility for SAC files
%
%     plot SAC seismograms read into matlab with rsac.m
% 
%     usage:  p1(file1, file2, ... , fileN, [time1 time2])
%
%     where file1,file2,... represents any matlab variable
%     that is a seismogram read in with rsac.m
%     [time1 time2] are optional.  If these are not specified
%     the entire time trace will be displayed, otherwise
%     all traces will be plotted over this time window.
%
%     The utility displays the following header variables:
%     
%     upper left:   'KSTNM'.'KCMPNM'
%     upper right:  'GCARC'
%     lower left:   name of matlab variable 
%
%     Markers are drawn if times are set for the 'O', and
%     'T0-T9' header variables.  These markers are also
%     labeled if 'KT0-KT9' are also set.  
%
%     Example:
%
%     To plot aak, casy0, casy1, and hrv for the entire time
%     series:
%
%     p1(aak,casy0,casy1,hrv) 
%
%     To just plot aak, and hrv for the time range of 0 to 300
%     seconds:
%
%     p1(aak,hrv,[0 300]) 
%
%     by Michael Thorne (5/2004)   mthorne@asu.edu

function p1(varargin) 

%set(gcf,'Name','P1 -- SAC Seismogram Plotting Utility', ...
%    'NumberTitle','off','Color',[.8 .8 .8], ...
%    'Pointer','crosshair','PaperOrientation','landscape', ...
%    'PaperPosition',[.5 .5 10 7.5],'PaperType','usletter');
set(gcf,'Name','P1 -- SAC Seismogram Plotting Utility', ...
    'NumberTitle','off','Color',[.8 .8 .8], ...
    'Pointer','crosshair', ...
    'PaperPosition',[.5 .5 10 7.5],'PaperType','usletter');

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


for i=1:loopend
  file=varargin{i};

  if nargin > 1 & i < loopend 
    subplot(loopend,1,i)
    plot(file(:,1),file(:,2))
    set(gca,'Xcolor',[.1 .1 .1],'Ycolor',[.1 .1 .1], ...
        'FontName','times','FontWeight','light', ...
        'TickDir','out','FontSize',8  )
    grid on
    axis tight
    if junk == 1
      ylimm=get(gca,'YLim');
      ymin=ylimm(1,1);
      ymax=ylimm(1,2);
      axis([xaxmin xaxmax ymin ymax])
    end

    %get axis limits
    ylimm=get(gca,'YLim');
    xlimm=get(gca,'XLim');
    ymin=ylimm(1,1); xmin=xlimm(1,1);
    ymax=ylimm(1,2); xmax=xlimm(1,2);
    yoff=.01*(ymax-ymin);
    xoff=.005*(xmax-xmin);

    hold on
    ot=lh(file,'O');
    if ot ~= -12345 & ot >= xmin & ot <= xmax
      plot([ot ot], [ymin ymax], 'Color',[1 0 0],'Linewidth',1)
      text(ot+xoff, ymin+yoff, 'O','Color',[1 0 0],'FontName','times', ...
       'VerticalAlignment','bottom','FontSize',8)
    end
    ot=lh(file,'T0');
    if ot ~= -12345 & ot >= xmin & ot <= xmax
      plot([ot ot], [ymin ymax], 'Color',[1 0 0],'Linewidth',1)
      phase=lh(file,'KT0');
      ans=strcmp(deblank(phase),'-12345');
      if (ans == 1)
       phase='T0';
      end
      text(ot+xoff, ymin+yoff, phase,'Color',[1 0 0],'FontName','times', ...
       'VerticalAlignment','bottom','FontSize',8)
    end
    ot=lh(file,'T1');
    if ot ~= -12345 & ot >= xmin & ot <= xmax
      plot([ot ot], [ymin ymax], 'Color',[1 0 0],'Linewidth',1)
      phase=lh(file,'KT1');
      ans=strcmp(deblank(phase),'-12345');
      if (ans == 1)
       phase='T1';
      end
      text(ot+xoff, ymin+yoff, phase,'Color',[1 0 0],'FontName','times', ...
       'VerticalAlignment','bottom','FontSize',8)
    end
    ot=lh(file,'T2');
    if ot ~= -12345 & ot >= xmin & ot <= xmax
      plot([ot ot], [ymin ymax], 'Color',[1 0 0],'Linewidth',1)
      phase=lh(file,'KT2');
      ans=strcmp(deblank(phase),'-12345');
      if (ans == 1)
       phase='T2';
      end
      text(ot+xoff, ymin+yoff, phase,'Color',[1 0 0],'FontName','times', ...
       'VerticalAlignment','bottom','FontSize',8)
    end
    ot=lh(file,'T3');
    if ot ~= -12345 & ot >= xmin & ot <= xmax
      plot([ot ot], [ymin ymax], 'Color',[1 0 0],'Linewidth',1)
      phase=lh(file,'KT3');
      ans=strcmp(deblank(phase),'-12345');
      if (ans == 1)
       phase='T3';
      end
      text(ot+xoff, ymin+yoff, phase,'Color',[1 0 0],'FontName','times', ...
       'VerticalAlignment','bottom','FontSize',8)
    end
    ot=lh(file,'T4');
    if ot ~= -12345 & ot >= xmin & ot <= xmax
      plot([ot ot], [ymin ymax], 'Color',[1 0 0],'Linewidth',1)
      phase=lh(file,'KT4');
      ans=strcmp(deblank(phase),'-12345');
      if (ans == 1)
       phase='T4';
      end
      text(ot+xoff, ymin+yoff, phase,'Color',[1 0 0],'FontName','times', ...
       'VerticalAlignment','bottom','FontSize',8)
    end
    ot=lh(file,'T5');
    if ot ~= -12345 & ot >= xmin & ot <= xmax
      plot([ot ot], [ymin ymax], 'Color',[1 0 0],'Linewidth',1)
      phase=lh(file,'KT5');
      ans=strcmp(deblank(phase),'-12345');
      if (ans == 1)
       phase='T5';
      end
      text(ot+xoff, ymin+yoff, phase,'Color',[1 0 0],'FontName','times', ...
       'VerticalAlignment','bottom','FontSize',8)
    end
    ot=lh(file,'T6');
    if ot ~= -12345 & ot >= xmin & ot <= xmax
      plot([ot ot], [ymin ymax], 'Color',[1 0 0],'Linewidth',1)
      phase=lh(file,'KT6');
      ans=strcmp(deblank(phase),'-12345');
      if (ans == 1)
       phase='T6';
      end
      text(ot+xoff, ymin+yoff, phase,'Color',[1 0 0],'FontName','times', ...
       'VerticalAlignment','bottom','FontSize',8)
    end
    ot=lh(file,'T7');
    if ot ~= -12345 & ot >= xmin & ot <= xmax
      plot([ot ot], [ymin ymax], 'Color',[1 0 0],'Linewidth',1)
      phase=lh(file,'KT7');
      ans=strcmp(deblank(phase),'-12345');
      if (ans == 1)
       phase='T7';
      end
      text(ot+xoff, ymin+yoff, phase,'Color',[1 0 0],'FontName','times', ...
       'VerticalAlignment','bottom','FontSize',8)
    end
    ot=lh(file,'T8');
    if ot ~= -12345 & ot >= xmin & ot <= xmax
      plot([ot ot], [ymin ymax], 'Color',[1 0 0],'Linewidth',1)
      phase=lh(file,'KT8');
      ans=strcmp(deblank(phase),'-12345');
      if (ans == 1)
       phase='T8';
      end
      text(ot+xoff, ymin+yoff, phase,'Color',[1 0 0],'FontName','times', ...
       'VerticalAlignment','bottom','FontSize',8)
    end
    ot=lh(file,'T9');
    if ot ~= -12345 & ot >= xmin & ot <= xmax
      plot([ot ot], [ymin ymax], 'Color',[1 0 0],'Linewidth',1)
      phase=lh(file,'KT9');
      ans=strcmp(deblank(phase),'-12345');
      if (ans == 1)
       phase='T9';
      end
      text(ot+xoff, ymin+yoff, phase,'Color',[1 0 0],'FontName','times', ...
       'VerticalAlignment','bottom','FontSize',8)
    end
    hold off

    t1=lh(file,'KSTNM');
    t3=lh(file,'KCMPNM');
    ans=strcmp(deblank(t3),'-12345');
    if (ans == 0)
     tname=strcat(t1,'.',t3);
    else
     tname=t1;
    end
    text(xmin+xoff, ymax-yoff, tname, 'Color',[0 .6 .9],'FontName','times', ...
       'VerticalAlignment','top','FontSize',8)
    t1=num2str(lh(file,'GCARC'));
    tdist=strcat(t1,'^o');
    % zw modified
    %text(xmax-xoff, ymax-yoff, tdist, 'Color',[0 .6 .9],'FontName','times', ...
    %   'VerticalAlignment','top','HorizontalAlignment','right', ...
    %   'FontSize',8)
    
    tname=inputname(i);
    text(xmin+xoff, ymin+yoff, tname, 'Color',[0 .6 .9],'FontName','times', ...
       'VerticalAlignment','bottom','HorizontalAlignment','left', ...
       'FontSize',8)
  
  else

    subplot(loopend,1,i)
    plot(file(:,1),file(:,2))
    set(gca,'Xcolor',[.1 .1 .1],'Ycolor',[.1 .1 .1], ...
        'FontName','times','FontWeight','light', ...
        'TickDir','out','FontSize',8  )
    grid on
    axis tight
    if junk == 1
      ylimm=get(gca,'YLim');
      ymin=ylimm(1,1);
      ymax=ylimm(1,2);
      axis([xaxmin xaxmax ymin ymax])
    end
    xlabel('Time (sec)')
  
    %get axis limits
    ylimm=get(gca,'YLim'); xlimm=get(gca,'XLim');
    ymin=ylimm(1,1); xmin=xlimm(1,1);
    ymax=ylimm(1,2); xmax=xlimm(1,2);
    yoff=.01*(ymax-ymin);
    xoff=.005*(xmax-xmin);

    hold on
    ot=lh(file,'O');
    if ot ~= -12345 & ot >= xmin & ot <= xmax
      plot([ot ot], [ymin ymax], 'Color',[1 0 0],'Linewidth',1)
      text(ot+xoff, ymin+yoff, 'O','Color',[1 0 0],'FontName','times', ...
       'VerticalAlignment','bottom','FontSize',8)
    end
    ot=lh(file,'T0');
    if ot ~= -12345 & ot >= xmin & ot <= xmax
      plot([ot ot], [ymin ymax], 'Color',[1 0 0],'Linewidth',1)
      phase=lh(file,'KT0');
      ans=strcmp(deblank(phase),'-12345');
      if (ans == 1)
       phase='T0';
      end
      text(ot+xoff, ymin+yoff, phase,'Color',[1 0 0],'FontName','times', ...
       'VerticalAlignment','bottom','FontSize',8)
    end
    ot=lh(file,'T1');
    if ot ~= -12345 & ot >= xmin & ot <= xmax
      plot([ot ot], [ymin ymax], 'Color',[1 0 0],'Linewidth',1)
      phase=lh(file,'KT1');
      ans=strcmp(deblank(phase),'-12345');
      if (ans == 1)
       phase='T1';
      end
      text(ot+xoff, ymin+yoff, phase,'Color',[1 0 0],'FontName','times', ...
       'VerticalAlignment','bottom','FontSize',8)
    end
    ot=lh(file,'T2');
    if ot ~= -12345 & ot >= xmin & ot <= xmax
      plot([ot ot], [ymin ymax], 'Color',[1 0 0],'Linewidth',1)
      phase=lh(file,'KT2');
      ans=strcmp(deblank(phase),'-12345');
      if (ans == 1)
       phase='T2';
      end
      text(ot+xoff, ymin+yoff, phase,'Color',[1 0 0],'FontName','times', ...
       'VerticalAlignment','bottom','FontSize',8)
    end
    ot=lh(file,'T3');
    if ot ~= -12345 & ot >= xmin & ot <= xmax
      plot([ot ot], [ymin ymax], 'Color',[1 0 0],'Linewidth',1)
      phase=lh(file,'KT3');
      ans=strcmp(deblank(phase),'-12345');
      if (ans == 1)
       phase='T3';
      end
      text(ot+xoff, ymin+yoff, phase,'Color',[1 0 0],'FontName','times', ...
       'VerticalAlignment','bottom','FontSize',8)
    end
    ot=lh(file,'T4');
    if ot ~= -12345 & ot >= xmin & ot <= xmax
      plot([ot ot], [ymin ymax], 'Color',[1 0 0],'Linewidth',1)
      phase=lh(file,'KT4');
      ans=strcmp(deblank(phase),'-12345');
      if (ans == 1)
       phase='T4';
      end
      text(ot+xoff, ymin+yoff, phase,'Color',[1 0 0],'FontName','times', ...
       'VerticalAlignment','bottom','FontSize',8)
    end
    ot=lh(file,'T5');
    if ot ~= -12345 & ot >= xmin & ot <= xmax
      plot([ot ot], [ymin ymax], 'Color',[1 0 0],'Linewidth',1)
      phase=lh(file,'KT5');
      ans=strcmp(deblank(phase),'-12345');
      if (ans == 1)
       phase='T5';
      end
      text(ot+xoff, ymin+yoff, phase,'Color',[1 0 0],'FontName','times', ...
       'VerticalAlignment','bottom','FontSize',8)
    end
    ot=lh(file,'T6');
    if ot ~= -12345 & ot >= xmin & ot <= xmax
      plot([ot ot], [ymin ymax], 'Color',[1 0 0],'Linewidth',1)
      phase=lh(file,'KT6');
      ans=strcmp(deblank(phase),'-12345');
      if (ans == 1)
       phase='T6';
      end
      text(ot+xoff, ymin+yoff, phase,'Color',[1 0 0],'FontName','times', ...
       'VerticalAlignment','bottom','FontSize',8)
    end
    ot=lh(file,'T7');
    if ot ~= -12345 & ot >= xmin & ot <= xmax
      plot([ot ot], [ymin ymax], 'Color',[1 0 0],'Linewidth',1)
      phase=lh(file,'KT7');
      ans=strcmp(deblank(phase),'-12345');
      if (ans == 1)
       phase='T7';
      end
      text(ot+xoff, ymin+yoff, phase,'Color',[1 0 0],'FontName','times', ...
       'VerticalAlignment','bottom','FontSize',8)
    end
    ot=lh(file,'T8');
    if ot ~= -12345 & ot >= xmin & ot <= xmax
      plot([ot ot], [ymin ymax], 'Color',[1 0 0],'Linewidth',1)
      phase=lh(file,'KT8');
      ans=strcmp(deblank(phase),'-12345');
      if (ans == 1)
       phase='T8';
      end
      text(ot+xoff, ymin+yoff, phase,'Color',[1 0 0],'FontName','times', ...
       'VerticalAlignment','bottom','FontSize',8)
    end
    ot=lh(file,'T9');
    if ot ~= -12345 & ot >= xmin & ot <= xmax
      plot([ot ot], [ymin ymax], 'Color',[1 0 0],'Linewidth',1)
      phase=lh(file,'KT9');
      ans=strcmp(deblank(phase),'-12345');
      if (ans == 1)
       phase='T9';
      end
      text(ot+xoff, ymin+yoff, phase,'Color',[1 0 0],'FontName','times', ...
       'VerticalAlignment','bottom','FontSize',8)
    end
    hold off

    t1=lh(file,'KSTNM');
    t3=lh(file,'KCMPNM');
    ans=strcmp(deblank(t3),'-12345');
    if (ans == 0)
     tname=strcat(t1,'.',t3);
    else
     tname=t1;
    end
    text(xmin+xoff, ymax-yoff, tname, 'Color',[0 .6 .9],'FontName','times', ...
       'VerticalAlignment','top','FontSize',8)
    t1=num2str(lh(file,'GCARC'));
    tdist=strcat(t1,'^o');
    % zw modified
    %text(xmax-xoff, ymax-yoff, tdist, 'Color',[0 .6 .9],'FontName','times', ...
    %   'VerticalAlignment','top','HorizontalAlignment','right', ...
    %   'FontSize',8)
    tname=inputname(i);
    text(xmin+xoff, ymin+yoff, tname, 'Color',[0 .6 .9],'FontName','times', ...
       'VerticalAlignment','bottom','HorizontalAlignment','left', ...
       'FontSize',8)

  end

end
