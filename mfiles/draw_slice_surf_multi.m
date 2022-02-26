% Draw the snapshot on multi slices by using surf style
% Author:   Yuanhang Huo
% Email:    yhhuo@mail.ustc.edu.cn
% Date:     2021.05.31

clear all;
addmypath
% -------------------------- parameters input -------------------------- %
% file and path name
parfnm='../project/test.json';
output_dir='../project/output';

% which slice to plot
% slice 1
slicedir{1}='y';
sliceid{1}=120;
% slice 2
slicedir{2}='x';
sliceid{2}=190;
% slice 3
slicedir{3}='z';
sliceid{3}=59;

% which variable and time to plot
varnm='Vz';
ns=1;
ne=500;
nt=50;

% figure control parameters
flag_km     = 1;
flag_emlast = 1;
flag_print  = 0;
scl_daspect =[1 1 1];
clrmp       = 'parula';
taut=0.5;
% ---------------------------------------------------------------------- %



% read parameter file
par=loadjson(parfnm);
ni=par.number_of_total_grid_points_x;
nj=par.number_of_total_grid_points_y;
nk=par.number_of_total_grid_points_z;
nproi=par.number_of_mpiprocs_x;
nproj=par.number_of_mpiprocs_y;


% figure plot
hid=figure;
set(hid,'BackingStore','on');

% load data
for nlayer=ns:nt:ne
    
    for i=1:length(sliceid)
        % -------------------- slice x ---------------------- %
        if slicedir{i} == 'x'
            
            for jp=0:nproj-1
                
                % snapshot data
                slicestruct{i}=dir([output_dir,'/','slicex_i',num2str(sliceid{i}),'_px*_py',num2str(jp),'.nc']);
                slicenm{i}=slicestruct{i}.name;
                slicenm{i}=[output_dir,'/',slicenm{i}];
                pnjstruct{i}=nc_getdiminfo(slicenm{i},'j');
                pnj{i}=pnjstruct{i}.Length;
                pnkstruct{i}=nc_getdiminfo(slicenm{i},'k');
                pnk{i}=pnkstruct{i}.Length;
                if jp==0
                    V{i}=squeeze(nc_varget(slicenm{i},varnm,[nlayer-1,0,0],[1,pnk{i},pnj{i}],[1,1,1]));
                else
                    V0{i}=squeeze(nc_varget(slicenm{i},varnm,[nlayer-1,0,0],[1,pnk{i},pnj{i}],[1,1,1]));
                    V{i}=horzcat(V{i},V0{i});
                end
                t=nc_varget(slicenm{i},'time',[nlayer-1],[1]);
                
                % coordinate data
                slicestruct{i}=dir([output_dir,'/','slicex_i',num2str(sliceid{i}),'_px*_py',num2str(jp),'.nc']);
                slicenm{i}=slicestruct{i}.name;
                ip=str2num(slicenm{i}( strfind(slicenm{i},'px')+2 : strfind(slicenm{i},'_py')-1 ));
                slicenm{i}=[output_dir,'/',slicenm{i}];
                coordnm{i}=['coord','_px',num2str(ip),'_py',num2str(jp),'.nc'];
                coordnm{i}=[output_dir,'/',coordnm{i}];
                idwithghost=double(nc_attget(slicenm{i},nc_global,'i_index_with_ghosts_in_this_thread'));
                coorddimstruct{i}=nc_getdiminfo(coordnm{i},'k');
                slicedimstruct{i}=nc_getdiminfo(slicenm{i},'k');
                ghostp=(coorddimstruct{i}.Length-slicedimstruct{i}.Length)/2;
                if jp==0
                    X{i}=squeeze(nc_varget(coordnm{i},'x',[ghostp,ghostp,idwithghost],[pnk{i},pnj{i},1],[1,1,1]));
                    Y{i}=squeeze(nc_varget(coordnm{i},'y',[ghostp,ghostp,idwithghost],[pnk{i},pnj{i},1],[1,1,1]));
                    Z{i}=squeeze(nc_varget(coordnm{i},'z',[ghostp,ghostp,idwithghost],[pnk{i},pnj{i},1],[1,1,1]));
                else
                    X0{i}=squeeze(nc_varget(coordnm{i},'x',[ghostp,ghostp,idwithghost],[pnk{i},pnj{i},1],[1,1,1]));
                    Y0{i}=squeeze(nc_varget(coordnm{i},'y',[ghostp,ghostp,idwithghost],[pnk{i},pnj{i},1],[1,1,1]));
                    Z0{i}=squeeze(nc_varget(coordnm{i},'z',[ghostp,ghostp,idwithghost],[pnk{i},pnj{i},1],[1,1,1]));
                    X{i}=horzcat(X{i},X0{i});
                    Y{i}=horzcat(Y{i},Y0{i});
                    Z{i}=horzcat(Z{i},Z0{i});
                end
                
            end
            
            % unit
            str_unit='m';
            if flag_km
                X{i}=X{i}/1e3;
                Y{i}=Y{i}/1e3;
                Z{i}=Z{i}/1e3;
                str_unit='km';
            end
            
            surf(X{i},Y{i},Z{i},V{i});
            xlabel(['X axis (',str_unit,')']);
            ylabel(['Y axis (',str_unit,')']);
            zlabel(['Z axis (',str_unit,')']);
%             view(45,15);
            
            % -------------------- slice y ---------------------- %
        elseif slicedir{i} == 'y'
            
            for ip=0:nproi-1
                
                % snapshot data
                slicestruct{i}=dir([output_dir,'/','slicey_j',num2str(sliceid{i}),'_px',num2str(ip),'_py*.nc']);
                slicenm{i}=slicestruct{i}.name;
                slicenm{i}=[output_dir,'/',slicenm{i}];
                pnistruct{i}=nc_getdiminfo(slicenm{i},'i');
                pni{i}=pnistruct{i}.Length;
                pnkstruct{i}=nc_getdiminfo(slicenm{i},'k');
                pnk{i}=pnkstruct{i}.Length;
                if ip==0
                    V{i}=squeeze(nc_varget(slicenm{i},varnm,[nlayer-1,0,0],[1,pnk{i},pni{i}],[1,1,1]));
                else
                    V0{i}=squeeze(nc_varget(slicenm{i},varnm,[nlayer-1,0,0],[1,pnk{i},pni{i}],[1,1,1]));
                    V{i}=horzcat(V{i},V0{i});
                end
                t=nc_varget(slicenm{i},'time',[nlayer-1],[1]);
                
                % coordinate data
                slicestruct{i}=dir([output_dir,'/','slicey_j',num2str(sliceid{i}),'_px',num2str(ip),'_py*.nc']);
                slicenm{i}=slicestruct{i}.name;
                jp=str2num(slicenm{i}( strfind(slicenm{i},'py')+2 : strfind(slicenm{i},'.nc')-1 ));
                slicenm{i}=[output_dir,'/',slicenm{i}];
                coordnm{i}=['coord','_px',num2str(ip),'_py',num2str(jp),'.nc'];
                coordnm{i}=[output_dir,'/',coordnm{i}];
                idwithghost=double(nc_attget(slicenm{i},nc_global,'j_index_with_ghosts_in_this_thread'));
                coorddimstruct{i}=nc_getdiminfo(coordnm{i},'i');
                slicedimstruct{i}=nc_getdiminfo(slicenm{i},'i');
                ghostp=(coorddimstruct{i}.Length-slicedimstruct{i}.Length)/2;
                if ip==0
                    X{i}=squeeze(nc_varget(coordnm{i},'x',[ghostp,idwithghost,ghostp],[pnk{i},1,pni{i}],[1,1,1]));
                    Y{i}=squeeze(nc_varget(coordnm{i},'y',[ghostp,idwithghost,ghostp],[pnk{i},1,pni{i}],[1,1,1]));
                    Z{i}=squeeze(nc_varget(coordnm{i},'z',[ghostp,idwithghost,ghostp],[pnk{i},1,pni{i}],[1,1,1]));
                else
                    X0{i}=squeeze(nc_varget(coordnm{i},'x',[ghostp,idwithghost,ghostp],[pnk{i},1,pni{i}],[1,1,1]));
                    Y0{i}=squeeze(nc_varget(coordnm{i},'y',[ghostp,idwithghost,ghostp],[pnk{i},1,pni{i}],[1,1,1]));
                    Z0{i}=squeeze(nc_varget(coordnm{i},'z',[ghostp,idwithghost,ghostp],[pnk{i},1,pni{i}],[1,1,1]));
                    X{i}=horzcat(X{i},X0{i});
                    Y{i}=horzcat(Y{i},Y0{i});
                    Z{i}=horzcat(Z{i},Z0{i});
                end
                
            end
            
            % unit
            str_unit='m';
            if flag_km
                X{i}=X{i}/1e3;
                Y{i}=Y{i}/1e3;
                Z{i}=Z{i}/1e3;
                str_unit='km';
            end
            
            surf(X{i},Y{i},Z{i},V{i});
            xlabel(['X axis (',str_unit,')']);
            ylabel(['Y axis (',str_unit,')']);
            zlabel(['Z axis (',str_unit,')']);
%             view(45,15);
            
            % -------------------- slice z ---------------------- %
        else
            
            for jp=0:nproj-1
                for ip=0:nproi-1
                    
                    % snapshot data
                    slicenm{i}=[output_dir,'/','slicez_k',num2str(sliceid{i}),'_px',num2str(ip),'_py',num2str(jp),'.nc'];
                    pnistruct{i}=nc_getdiminfo(slicenm{i},'i');
                    pni{i}=pnistruct{i}.Length;
                    pnjstruct{i}=nc_getdiminfo(slicenm{i},'j');
                    pnj{i}=pnjstruct{i}.Length;
                    if ip==0
                        VV{i}=squeeze(nc_varget(slicenm{i},varnm,[nlayer-1,0,0],[1,pnj{i},pni{i}],[1,1,1]));
                    else
                        VV0{i}=squeeze(nc_varget(slicenm{i},varnm,[nlayer-1,0,0],[1,pnj{i},pni{i}],[1,1,1]));
                        VV{i}=horzcat(VV{i},VV0{i});
                    end
                    t=nc_varget(slicenm{i},'time',[nlayer-1],[1]);
                    
                    % coordinate data
                    coordnm{i}=['coord','_px',num2str(ip),'_py',num2str(jp),'.nc'];
                    coordnm{i}=[output_dir,'/',coordnm{i}];
                    idwithghost=double(nc_attget(slicenm{i},nc_global,'k_index_with_ghosts_in_this_thread'));
                    coorddimstruct{i}=nc_getdiminfo(coordnm{i},'j');
                    slicedimstruct{i}=nc_getdiminfo(slicenm{i},'j');
                    ghostp=(coorddimstruct{i}.Length-slicedimstruct{i}.Length)/2;
                    if ip==0
                        XX{i}=squeeze(nc_varget(coordnm{i},'x',[idwithghost,ghostp,ghostp],[1,pnj{i},pni{i}],[1,1,1]));
                        YY{i}=squeeze(nc_varget(coordnm{i},'y',[idwithghost,ghostp,ghostp],[1,pnj{i},pni{i}],[1,1,1]));
                        ZZ{i}=squeeze(nc_varget(coordnm{i},'z',[idwithghost,ghostp,ghostp],[1,pnj{i},pni{i}],[1,1,1]));
                    else
                        XX0{i}=squeeze(nc_varget(coordnm{i},'x',[idwithghost,ghostp,ghostp],[1,pnj{i},pni{i}],[1,1,1]));
                        YY0{i}=squeeze(nc_varget(coordnm{i},'y',[idwithghost,ghostp,ghostp],[1,pnj{i},pni{i}],[1,1,1]));
                        ZZ0{i}=squeeze(nc_varget(coordnm{i},'z',[idwithghost,ghostp,ghostp],[1,pnj{i},pni{i}],[1,1,1]));
                        XX{i}=horzcat(XX{i},XX0{i});
                        YY{i}=horzcat(YY{i},YY0{i});
                        ZZ{i}=horzcat(ZZ{i},ZZ0{i});
                    end
                    
                end
                
                if jp==0
                    V{i}=VV{i};
                    X{i}=XX{i};
                    Y{i}=YY{i};
                    Z{i}=ZZ{i};
                else
                    V{i}=vertcat(V{i},VV{i});
                    X{i}=vertcat(X{i},XX{i});
                    Y{i}=vertcat(Y{i},YY{i});
                    Z{i}=vertcat(Z{i},ZZ{i});
                end
            end
            
            % unit
            str_unit='m';
            if flag_km
                X{i}=X{i}/1e3;
                Y{i}=Y{i}/1e3;
                Z{i}=Z{i}/1e3;
                str_unit='km';
            end
            
            surf(X{i},Y{i},Z{i},V{i});
            xlabel(['X axis (',str_unit,')']);
            ylabel(['Y axis (',str_unit,')']);
            zlabel(['Z axis (',str_unit,')']);
            
        end
    
        hold on;
        
    end
    
    hold off;
    
    disp([ '  draw ' num2str(nlayer) 'th time step (t=' num2str(t) ')']);
    
    set(gca,'layer','top');
    set(gcf,'color','white','renderer','painters');

    % axis image
    % shading interp;
    shading flat;
    % colorbar range/scale
    if exist('scl_caxis')
        caxis(scl_caxis);
    end
    % axis daspect
    if exist('scl_daspect')
        daspect(scl_daspect);
    end
    % colormap and colorbar
    if exist('clrmp')
        colormap(clrmp);
    end
    colorbar('vert');
    
    % title
    titlestr=['Snapshot of ' varnm ' at ' ...
              '{\fontsize{12}{\bf ' ...
              num2str((t),'%7.3f') ...
              '}}s'];
    title(titlestr);
    
    drawnow;
    pause(taut);
    
    % save and print figure
    if flag_print==1
        width= 500;
        height=500;
        set(gcf,'paperpositionmode','manual');
        set(gcf,'paperunits','points');
        set(gcf,'papersize',[width,height]);
        set(gcf,'paperposition',[0,0,width,height]);
        fnm_out=[varnm '_ndim_',num2str(nlayer,'%5.5i')];
        print(gcf,[fnm_out '.png'],'-dpng');
    end
    
end
        
        
        