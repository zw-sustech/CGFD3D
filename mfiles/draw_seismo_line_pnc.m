% Plot the seismograms on a line

clear all;
%close all;
clc;
addmypath;

% -------------------------- parameters input -------------------------- %

%nc_file = '~/work/cgfd3d/line2/output/line_x_1.nc'
nc_file = '~/work/cgfd3d/line2/output/line_y_1.nc'

% which variable to plot
varnm='Vy';

% which station to plot (start from index '1')
ir_start = 0;
ir_end   = 20;
ir_skip  = 1;

%- derive
ir_count = (ir_end - ir_start) / ir_skip;

% figure control parameters
flag_print=0;

% ---------------------------------------------------------------------- %

% read var
T       = ncread(nc_file, 'time');
V       = ncread(nc_file, varnm, [1,ir_start+1], [inf,ir_count],[1,ir_skip]);

% plot receiver line

scl = max(max(abs(V)));

ytickincre=10;

figure;

for irec = 0 : ir_count - 1
    plot(T, V(:,irec+1) + irec*(2*scl), 'b');
    hold on;
end

xlabel('Time (s)');
ylabel('Receiver ID');
%title([varnm,' at No.',num2str(lineid), ' Receiver Line ','(',linenm,')'],'interpreter','none');
set(gca,'ytick',[0:10:ir_count-1]*(2*scl),'yticklabel',[1:ir_count]);
set(gca,'ylim',[-2,ir_count+1]*(2*scl));
set(gcf,'position',[0,0,400,1200]);
set(gcf,'color','white','renderer','painters');

% save and print figure
if flag_print
    width= 400;
    height=1200;
    set(gcf,'paperpositionmode','manual');
    set(gcf,'paperunits','points');
    set(gcf,'papersize',[width,height]);
    set(gcf,'paperposition',[0,0,width,height]);
    %print(gcf,[varnm,'_line_no',num2str(lineid),'.png'],'-dpng');
end

