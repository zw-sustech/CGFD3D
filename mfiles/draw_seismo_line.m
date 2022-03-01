% Plot the seismograms on a line
% Author:   Yuanhang Huo
% Email:    yhhuo@mail.ustc.edu.cn
% Date:     2021.06.04

clear all;
addmypath;
% -------------------------- parameters input -------------------------- %
% file and path name
parfnm='../project/test.json';
output_dir='../project/output';

% which line to plot (start from index '1')
lineid=1;

% which receiver of the line to plot (start from index '0')
recid=1;

% which variable to plot
varnm='Vz';

% figure control parameters
flag_print=1;

% ---------------------------------------------------------------------- %

% read parameter file
par=loadjson(parfnm);

fileID = fopen(par.in_source_file);
lineprefix = fgetl(fileID);
while(lineprefix(1) == "#")
    lineprefix = fgetl(fileID);
    if(lineprefix(1) ~= "#")

        break;
    end
end

% line name and receiver number
linenm=par.receiver_line{lineid}.name;
nrec=par.receiver_line{lineid}.grid_index_count;

% load data
for irec=0:nrec-1
    
    sacnm=[output_dir,'/',lineprefix,'.',linenm,'.','no',num2str(irec),'.',varnm,'.sac'];
    sacdata=rsac(sacnm);
    seismodata(irec+1,:)=sacdata(:,2);
    seismot(irec+1,:)=sacdata(:,1);
    
end


% plot single receiver
figure;
plot(seismot(recid+1,:),seismodata(recid+1,:),'b','linewidth',1.0);
xlabel('Time (s)');
ylabel('Amplitude');
title([varnm, ' at No.',num2str(recid+1),' Receiver of No.',num2str(lineid),...
       ' Line ','(',linenm,')'],'interpreter','none');
set(gcf,'color','white','renderer','painters');

% save and print figure
if flag_print
    width= 800;
    height=400;
    set(gcf,'paperpositionmode','manual');
    set(gcf,'paperunits','points');
    set(gcf,'papersize',[width,height]);
    set(gcf,'paperposition',[0,0,width,height]);
    print(gcf,[varnm,'_rec_no',num2str(recid),'.png'],'-dpng');
end


% plot receiver line
scl=max(max(abs(seismodata)));
% increment of show for Y-axis ticklabel
ytickincre=10;
figure;
for irec=0:nrec-1
    plot(seismot(irec+1,:),seismodata(irec+1,:)+irec*(2*scl),'b');
    hold on;
end
xlabel('Time (s)');
ylabel('Receiver ID');
title([varnm,' at No.',num2str(lineid), ' Receiver Line ','(',linenm,')'],'interpreter','none');
set(gca,'ytick',[0:10:nrec-1]*(2*scl),'yticklabel',[1:10:nrec]);
set(gca,'ylim',[-2,nrec+1]*(2*scl));
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
    print(gcf,[varnm,'_line_no',num2str(lineid),'.png'],'-dpng');
end
    

