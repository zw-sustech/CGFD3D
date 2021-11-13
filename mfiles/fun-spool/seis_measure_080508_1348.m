function [varargout]=seis_measure(src,event,varargin)
% help should be here

%  Initial variable and flag
V=varargin; nVar=nargin-2;

titlestr=[]; 
for n=1:nVar
    V{n}.name=strtrim(V{n}.name);
    titlestr=[titlestr ' ' V{n}.name];
    %liststr{n}=V{n}.spec_line;
    liststr{n}=V{n}.name;
    Vpos(n)=0; Vamp(n)=1; Vtime(n)=0;
end

%  Initialization tasks
fid=figure('position',[700 100 800 600], 'renderer','zbuffer', ...
     'toolbar','figure','menubar','figure');
set(fid,'Name',[titlestr ' for ' V{1}.KEVNM]);

set(fid,'visible','off')

ha = axes('Units','normalized','Position',[0.05,0.11,0.83,0.815]);

if 1
for n=1:nVar
    hold on
    pid(n)=plot(V{n}.t,V{n}.v/V{n}.v0,V{n}.spec_line);
end

scl_ylim=get(gca,'ylim'); scl_xlim=get(gca,'xlim');
Lx=(scl_xlim(2)-scl_xlim(1))/100; Ly=(scl_ylim(2)-scl_ylim(1))/11;

set(gca,'ylim',[scl_ylim(1)-Ly,scl_ylim(2)+Ly]);

y0=scl_ylim(2); y1=scl_ylim(2)-Ly;
x1=2*Lx; x2=x1+2*Lx; x0=x2+Lx;

for n=1:nVar
    plot([x1 x2],[y0 y0]-(n-1)*Ly,V{n}.spec_line);
    text(x0,y0-(n-1)*Ly,[V{n}.name '(' num2str(V{n}.v0,'%3.2e') ')']);
end

if isfield(V{1},'title'), title(V{1}.title); end

%  Construct the components.
if 1
%bgh_var = uibuttongroup('Parent',fid,'Title','Variable',...
%            'Position',[.9 .02 .08 .06]);
%rbh_var_R = uicontrol(bgh_var,'Style','radiobutton','String','R',...
%                'Units','normalized',...
%                'Position',[.05 .1 .45 .8]);
%rbh_var_S = uicontrol(bgh_var,'Style','radiobutton','String','S',...
%                'Units','normalized',...
%                'Position',[.5 .1 .45 .8]);

lbh = uicontrol(fid,'Style','popupmenu',...
                'String',liststr,...
                'fontsize',8, ...
                'Value',1,'Unit','normalized','Position',[.9 .02 .09 .032], ...
                'Callback',{@set_slider_pos});

bgh_tsk = uibuttongroup('Parent',fid,'Title','Move',...
            'Position',[.9 .10 .08 .12], ...
            'SelectionChangeFcn',{@set_slider_pos});
rbh_tsk_a = uicontrol(bgh_tsk,'Style','radiobutton','String','amplitude',...
                'Units','normalized',...
                'Position',[.1 .6 .8 .2]);
rbh_tsk_p = uicontrol(bgh_tsk,'Style','radiobutton','String','position',...
                'Units','normalized',...
                'Position',[.1 .2 .8 .2]);

hslider_x = uicontrol(fid,'Style','slider',...
                      'Max',50,'Min',0,'Value',0, ...
                      'SliderStep',[0.005 0.1], ...
                      'Unit','normalized','Position',[0.05 0.02 0.83 0.03], ...
                      'Callback',{@shift_time_axis});

hslider_y = uicontrol(fid,'Style','slider',...
                      'Max',1,'Min',-1,'Value',0, ...
                      'SliderStep',[0.01 0.1], ...
                      'Unit','normalized','Position',[0.92 0.25 0.02 0.68], ...
                      'Callback',{@modify_amplitude});
end

set(fid,'visible','on')

if nargout>0, varargout{1}=pid(1); end
if nargout>1, varargout{2}=pid(2); end
end

%  Callbacks
% ---------------------------------------------------------------------------

function shift_time_axis(hObject, eventdata)
v0= get(hObject,'Value');
n=get(lbh,'Value');
set(pid(n),'xdata',get(pid(n),'xdata')-Vtime(n)+v0);
Vtime(n)=v0;
end

% ---------------------------------------------------------------------------
function modify_amplitude(hObject, eventdata)
v0= get(hObject,'Value');
v1=get(hObject,'min'); v2=get(hObject,'max');
Ly=(v2-v1)/2/5;

n=get(lbh,'Value');

if ( get(rbh_tsk_a,'Value')==get(rbh_tsk_a,'Max') )
   if v0>=0, fct=1+v0/Ly; else, fct=1+v0; end
   set(pid(n),'ydata',(get(pid(n),'ydata')-Vpos(n))/Vamp(n)*fct+Vpos(n));
   Vamp(n)=fct;
else
   set(pid(n),'ydata',(get(pid(n),'ydata')-Vpos(n)+v0));
   Vpos(n)=v0;
end

end

% ---------------------------------------------------------------------------
function set_slider_pos(src, eventdata)
n=get(lbh,'Value');
% y axis slider
if (get(rbh_tsk_a,'Value')==get(rbh_tsk_a,'Max') )
   v1=get(hslider_y,'min'); v2=get(hslider_y,'max'); Ly=(v2-v1)/2/5;
   if Vamp(n)>=1, y0=(Vamp(n)-1)*Ly; else, y0=Vamp(n)-1; end
   set(hslider_y,'Value',y0);
else
   set(hslider_y,'Value',Vpos(n));
end
% x axis slider
set(hslider_x,'Value',Vtime(n));

end

% ---------------------------------------------------------------------------

end
