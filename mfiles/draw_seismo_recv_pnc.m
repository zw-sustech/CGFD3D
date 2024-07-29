% Plot the seismograms on a line

clear all;
%close all;
clc;
addmypath;

% -------------------------- parameters input -------------------------- %

nc_file = '~/work/cgfd3d/sta/output/station.nc'

% which variable to plot
varnm='Vx';

% which station to plot (start from index '1')
seq_no=5;

% figure control parameters
flag_print=0;

% ---------------------------------------------------------------------- %

% read var
T       = ncread(nc_file, 'time');
V       = ncread(nc_file, varnm, [1,seq_no], [inf,1]);
STANM_NC = ncread(nc_file, 'station_name', [1,seq_no], [inf,1]);

STANM = deblank(STANM_NC')

% plot receiver
for irec = seq_no : 1 : seq_no

    figure;

    plot(T,V,'b','linewidth',1.0);

    xlabel('Time (s)');
    ylabel('Amplitude');
    %title([varnm, ' recv No.',num2str(seq_no),' interpreter ','yes']);
    title([varnm, ' at recv ',STANM,' interpreter ','yes']);
    set(gcf,'color','white','renderer','painters');
    % save and print figure
    if flag_print
        width= 800;
        height=400;
        set(gcf,'paperpositionmode','manual');
        set(gcf,'paperunits','points');
        set(gcf,'papersize',[width,height]);
        set(gcf,'paperposition',[0,0,width,height]);
        print(gcf,[varnm,'_rec_no',num2str(irec),'.png'],'-dpng');
    end
end
