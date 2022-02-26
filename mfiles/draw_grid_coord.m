% Plot the grid
% Author:       Yuanhang Huo
% Email:        yhhuo@mail.ustc.edu.cn
% Affiliation:  University of Science and Technology of China
% Date:         2021.06.06

clear all;
addmypath
% -------------------------- parameters input -------------------------- %
% file and path name
parfnm='../project/test.json';
output_dir='../project/output';

% which grid profile to plot
subs=[5,5,1];     % start from index '1'
subc=[-1,-1,1];     % '-1' to plot all points in this dimension
subt=[1,1,1];

% figure control parameters
flag_km     = 1;
flag_emlast = 1;
flag_print  = 0;
flag_title  = 1;
scl_daspect = [1 1 1];
% ---------------------------------------------------------------------- %



% load grid coordinate
coordinfo=locate_coord(parfnm,'start',subs,'count',subc,'stride',subt,'coorddir',output_dir);
[x,y,z]=gather_coord(coordinfo,'coorddir',output_dir);
nx=size(x,1);
ny=size(x,2);
nz=size(x,3);
% coordinate unit
str_unit='m';
if flag_km
   x=x/1e3;
   y=y/1e3;
   z=z/1e3;
   str_unit='km';
end

% plot
x=squeeze(x);
y=squeeze(y);
z=squeeze(z);
pltincre1=2;
pltincre2=2;
figure;
if nx == 1
    plot(permute(y(1:pltincre1:end,1:pltincre2:end),[2,1]),...
         permute(z(1:pltincre1:end,1:pltincre2:end),[2,1]),...
        'k-');
    hold on;
    plot(y(1:pltincre1:end,1:pltincre2:end),...
         z(1:pltincre1:end,1:pltincre2:end),...
        'k-');
    xlabel(['Y axis (' str_unit ')']);
    ylabel(['Z axis (' str_unit ')']);
     
elseif ny == 1
    plot(permute(x(1:pltincre1:end,1:pltincre2:end),[2,1]),...
         permute(z(1:pltincre1:end,1:pltincre2:end),[2,1]),...
         'k-');
    hold on;
    plot(x(1:pltincre1:end,1:pltincre2:end),...
         z(1:pltincre1:end,1:pltincre2:end),...
         'k-');
    xlabel(['X axis (' str_unit ')']);
    ylabel(['Z axis (' str_unit ')']);
     
else
    plot(permute(x(1:pltincre1:end,1:pltincre2:end),[2,1]),...
         permute(y(1:pltincre1:end,1:pltincre2:end),[2,1]),...
         'k-');
    hold on;
    plot(x(1:pltincre1:end,1:pltincre2:end),...
         y(1:pltincre1:end,1:pltincre2:end),...
         'k-');
    xlabel(['X axis (' str_unit ')']);
    ylabel(['Y axis (' str_unit ')']);
end

set(gca,'layer','top');
set(gcf,'color','white','renderer','painters');

% axis daspect
if exist('scl_daspect')
    daspect(scl_daspect);
end
axis tight;

% title
if flag_title
    if nx == 1
        gridtitle='YOZ-Grid';
    elseif ny == 1
        gridtitle='XOZ-Grid';
    else
        gridtitle='XOY-Grid';
    end
    title(gridtitle);
end

% save and print figure
if flag_print
    width= 500;
    height=500;
    set(gcf,'paperpositionmode','manual');
    set(gcf,'paperunits','points');
    set(gcf,'papersize',[width,height]);
    set(gcf,'paperposition',[0,0,width,height]);
    print(gcf,[gridtitle '.png'],'-dpng');
end


