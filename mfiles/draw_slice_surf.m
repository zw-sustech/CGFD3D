% Draw the snapshot on a slice by using surf style
% Author:   Yuanhang Huo
% Email:    yhhuo@mail.ustc.edu.cn
% Date:     2021.05.31

clear all;

% -------------------------- parameters input -------------------------- %
% file and path name
parfnm='./project/test.json';
output_dir='./project/output';

% which slice to plot
slicedir='x';
sliceid=19;

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
set(hid,'renderer','painters');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf,'PaperUnits','points');
set(gcf,'PaperPosition',[0 0 800 800]);

% load data
for nlayer=ns:nt:ne

    % -------------------- slice x ---------------------- %
    if slicedir == 'x'
        
        for jp=0:nproj-1
            
            % snapshot data
            slicestruct=dir([output_dir,'/','slicex_i',num2str(sliceid),'_px*_py',num2str(jp),'.nc']);
            slicenm=slicestruct.name;
            slicenm=[output_dir,'/',slicenm];
            pnjstruct=nc_getdiminfo(slicenm,'j');
            pnj=pnjstruct.Length;
            pnkstruct=nc_getdiminfo(slicenm,'k');
            pnk=pnkstruct.Length;
            if jp==0
                V=squeeze(nc_varget(slicenm,varnm,[nlayer-1,0,0],[1,pnk,pnj],[1,1,1]));
            else
                V0=squeeze(nc_varget(slicenm,varnm,[nlayer-1,0,0],[1,pnk,pnj],[1,1,1]));
                V=horzcat(V,V0);
            end
            t=nc_varget(slicenm,'time',[nlayer-1],[1]);
            
            % coordinate data
            slicestruct=dir([output_dir,'/','slicex_i',num2str(sliceid),'_px*_py',num2str(jp),'.nc']);
            slicenm=slicestruct.name;
            ip=str2num(slicenm( strfind(slicenm,'px')+2 : strfind(slicenm,'_py')-1 ));
            slicenm=[output_dir,'/',slicenm];
            coordnm=['coord','_px',num2str(ip),'_py',num2str(jp),'.nc'];
            coordnm=[output_dir,'/',coordnm];
            idwithghost=double(nc_attget(slicenm,nc_global,'i_index_with_ghosts_in_this_thread'));
            coorddimstruct=nc_getdiminfo(coordnm,'k');
            slicedimstruct=nc_getdiminfo(slicenm,'k');
            ghostp=(coorddimstruct.Length-slicedimstruct.Length)/2;
            if jp==0
                X=squeeze(nc_varget(coordnm,'x',[ghostp,ghostp,idwithghost],[pnk,pnj,1],[1,1,1]));
                Y=squeeze(nc_varget(coordnm,'y',[ghostp,ghostp,idwithghost],[pnk,pnj,1],[1,1,1]));
                Z=squeeze(nc_varget(coordnm,'z',[ghostp,ghostp,idwithghost],[pnk,pnj,1],[1,1,1]));
            else
                X0=squeeze(nc_varget(coordnm,'x',[ghostp,ghostp,idwithghost],[pnk,pnj,1],[1,1,1]));
                Y0=squeeze(nc_varget(coordnm,'y',[ghostp,ghostp,idwithghost],[pnk,pnj,1],[1,1,1]));
                Z0=squeeze(nc_varget(coordnm,'z',[ghostp,ghostp,idwithghost],[pnk,pnj,1],[1,1,1]));
                X=horzcat(X,X0);
                Y=horzcat(Y,Y0);
                Z=horzcat(Z,Z0);
            end  
            
        end
        
        str_unit='m';
        if flag_km
            X=X/1e3;
            Y=Y/1e3;
            Z=Z/1e3;
            str_unit='km';
        end
            
        surf(X,Y,Z,V);
        xlabel(['X axis (',str_unit,')']);
        ylabel(['Y axis (',str_unit,')']);
        zlabel(['Z axis (',str_unit,')']);
        view(45,15);
    
    % -------------------- slice y ---------------------- %
    elseif slicedir == 'y'
   
        for ip=0:nproi-1
            
            % snapshot data
            slicestruct=dir([output_dir,'/','slicey_j',num2str(sliceid),'_px',num2str(ip),'_py*.nc']);
            slicenm=slicestruct.name;
            slicenm=[output_dir,'/',slicenm];
            pnistruct=nc_getdiminfo(slicenm,'i');
            pni=pnistruct.Length;
            pnkstruct=nc_getdiminfo(slicenm,'k');
            pnk=pnkstruct.Length;
            if ip==0
                V=squeeze(nc_varget(slicenm,varnm,[nlayer-1,0,0],[1,pnk,pni],[1,1,1]));
            else
                V0=squeeze(nc_varget(slicenm,varnm,[nlayer-1,0,0],[1,pnk,pni],[1,1,1]));
                V=horzcat(V,V0);
            end
            t=nc_varget(slicenm,'time',[nlayer-1],[1]);
            
            % coordinate data
            slicestruct=dir([output_dir,'/','slicey_j',num2str(sliceid),'_px',num2str(ip),'_py*.nc']);
            slicenm=slicestruct.name;
            jp=str2num(slicenm( strfind(slicenm,'py')+2 : strfind(slicenm,'.nc')-1 ));
            slicenm=[output_dir,'/',slicenm];
            coordnm=['coord','_px',num2str(ip),'_py',num2str(jp),'.nc'];
            coordnm=[output_dir,'/',coordnm];
            idwithghost=double(nc_attget(slicenm,nc_global,'j_index_with_ghosts_in_this_thread'));
            coorddimstruct=nc_getdiminfo(coordnm,'i');
            slicedimstruct=nc_getdiminfo(slicenm,'i');
            ghostp=(coorddimstruct.Length-slicedimstruct.Length)/2;
            if ip==0
                X=squeeze(nc_varget(coordnm,'x',[ghostp,idwithghost,ghostp],[pnk,1,pni],[1,1,1]));
                Y=squeeze(nc_varget(coordnm,'y',[ghostp,idwithghost,ghostp],[pnk,1,pni],[1,1,1]));
                Z=squeeze(nc_varget(coordnm,'z',[ghostp,idwithghost,ghostp],[pnk,1,pni],[1,1,1]));
            else
                X0=squeeze(nc_varget(coordnm,'x',[ghostp,idwithghost,ghostp],[pnk,1,pni],[1,1,1]));
                Y0=squeeze(nc_varget(coordnm,'y',[ghostp,idwithghost,ghostp],[pnk,1,pni],[1,1,1]));
                Z0=squeeze(nc_varget(coordnm,'z',[ghostp,idwithghost,ghostp],[pnk,1,pni],[1,1,1]));
                X=horzcat(X,X0);
                Y=horzcat(Y,Y0);
                Z=horzcat(Z,Z0);
            end  
            
        end
        
        str_unit='m';
        if flag_km
            X=X/1e3;
            Y=Y/1e3;
            Z=Z/1e3;
            str_unit='km';
        end
            
        surf(X,Y,Z,V);
        xlabel(['X axis (',str_unit,')']);
        ylabel(['Y axis (',str_unit,')']);
        zlabel(['Z axis (',str_unit,')']);
        view(45,15);
    
    % -------------------- slice z ---------------------- %
    else

        for jp=0:nproj-1
            for ip=0:nproi-1
                
                % snapshot data
                slicenm=[output_dir,'/','slicez_k',num2str(sliceid),'_px',num2str(ip),'_py',num2str(jp),'.nc'];
                pnistruct=nc_getdiminfo(slicenm,'i');
                pni=pnistruct.Length;
                pnjstruct=nc_getdiminfo(slicenm,'j');
                pnj=pnjstruct.Length;
                if ip==0
                    VV=squeeze(nc_varget(slicenm,varnm,[nlayer-1,0,0],[1,pnj,pni],[1,1,1]));
                else
                    VV0=squeeze(nc_varget(slicenm,varnm,[nlayer-1,0,0],[1,pnj,pni],[1,1,1]));
                    VV=horzcat(VV,VV0);
                end
                t=nc_varget(slicenm,'time',[nlayer-1],[1]);
                
                % coordinate data
                coordnm=['coord','_px',num2str(ip),'_py',num2str(jp),'.nc'];
                coordnm=[output_dir,'/',coordnm];
                idwithghost=double(nc_attget(slicenm,nc_global,'k_index_with_ghosts_in_this_thread'));
                coorddimstruct=nc_getdiminfo(coordnm,'j');
                slicedimstruct=nc_getdiminfo(slicenm,'j');
                ghostp=(coorddimstruct.Length-slicedimstruct.Length)/2;
                if ip==0
                    XX=squeeze(nc_varget(coordnm,'x',[idwithghost,ghostp,ghostp],[1,pnj,pni],[1,1,1]));
                    YY=squeeze(nc_varget(coordnm,'y',[idwithghost,ghostp,ghostp],[1,pnj,pni],[1,1,1]));
                    ZZ=squeeze(nc_varget(coordnm,'z',[idwithghost,ghostp,ghostp],[1,pnj,pni],[1,1,1]));
                else
                    XX0=squeeze(nc_varget(coordnm,'x',[idwithghost,ghostp,ghostp],[1,pnj,pni],[1,1,1]));
                    YY0=squeeze(nc_varget(coordnm,'y',[idwithghost,ghostp,ghostp],[1,pnj,pni],[1,1,1]));
                    ZZ0=squeeze(nc_varget(coordnm,'z',[idwithghost,ghostp,ghostp],[1,pnj,pni],[1,1,1]));
                    XX=horzcat(XX,XX0);
                    YY=horzcat(YY,YY0);
                    ZZ=horzcat(ZZ,ZZ0);
                end
                    
            end
            
            if jp==0
                V=VV;
                X=XX;
                Y=YY;
                Z=ZZ;
            else
                V=vertcat(V,VV);
                X=vertcat(X,XX);
                Y=vertcat(Y,YY);
                Z=vertcat(Z,ZZ);
            end
        end
        
        str_unit='m';
        if flag_km
            X=X/1e3;
            Y=Y/1e3;
            Z=Z/1e3;
            str_unit='km';
        end
            
        surf(X,Y,Z,V);
        xlabel(['X axis (',str_unit,')']);
        ylabel(['Y axis (',str_unit,')']);
        zlabel(['Z axis (',str_unit,')']);
        
    end
    
    disp([ '  draw ' num2str(nlayer) 'th time step (t=' num2str(t) ')']);
    
    set(gca,'layer','top');

    %axis image
    %shading interp;
    shading flat;
    if exist('scl_caxis')
        caxis(scl_caxis);
    end
    if exist('scl_daspect')
        daspect(scl_daspect);
    end
    if exist('clrmp')
        colormap(clrmp);
    end
    colorbar('vert');
    
    titlestr=['Snapshot of ' varnm ' at ' ...
              '{\fontsize{12}{\bf ' ...
              num2str((t),'%7.3f') ...
              '}}s'];
    title(titlestr);
    
    drawnow;
    pause(taut);
    
    if flag_print==1
        fnm_out=[varnm '_ndim',num2str(nlayer,'%5.5i')];
%         set(gca,'FontName','FixedWidth');
        print(gcf,[fnm_out '.png'],'-dpng');
    end
    
end
        
        
        
        
        

        
