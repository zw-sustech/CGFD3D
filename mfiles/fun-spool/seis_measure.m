function [varargout]=seis_measure(src,event,varargin)
% help should be here

%  Initial variable and flag
flag_phase=0;
nVar=nargin-2;
%V=varargin;
V{1}=varargin{1};
idx=1;
for n=2:nVar
    for m=1:numel(varargin{n})
        idx=idx+1;
        V{idx}=varargin{n}{m};
    end
end
nVar=idx;

for n=1:nVar
    if isfield(V{n},'phase')
       flag_phase=1; P=V{n}.phase; break;
    end
end

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

[hlegd_line,hlegd_text]=plot_legend;
[phspid,phstid]=plot_phase;

if isfield(V{1},'title'), title(V{1}.title); end

%  Construct the components.
% -- popupmenu --
lbh = uicontrol(fid,'Style','popupmenu',...
                'String',liststr,...
                'fontsize',8, ...
                'Value',1,'Unit','normalized','Position',[.9 .02 .09 .032], ...
                'Callback',{@set_slider_pos});

% -- ratio button --
bgh_tsk = uibuttongroup('Parent',fid,'Title','Move',...
            'Position',[.9 .10 .08 .12], ...
            'SelectionChangeFcn',{@set_slider_pos});
rbh_tsk_a = uicontrol(bgh_tsk,'Style','radiobutton','String','amplitude',...
                'Units','normalized',...
                'Position',[.1 .6 .8 .2]);
rbh_tsk_p = uicontrol(bgh_tsk,'Style','radiobutton','String','position',...
                'Units','normalized',...
                'Position',[.1 .2 .8 .2]);

% -- slider --
htext_x1 = uicontrol(fid,'Style','text', ...
           'String','0','Unit','normalized','Position',[0.05 0.001 0.02 0.02], ...
           'backgroundcolor',get(fid,'color'));
htext_x2 = uicontrol(fid,'Style','text', ...
           'String','50','Unit','normalized','Position',[0.05+0.83-0.02 0.001 0.02 0.02], ...
           'backgroundcolor',get(fid,'color'));
htext_x = uicontrol(fid,'Style','text', ...
           'String','0','Unit','normalized','Position',[0.05 0.001 0.02 0.02], ...
           'backgroundcolor',get(fid,'color'), ...
           'visible','off');
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
% -- zoom --
hzoom=zoom(gca);
set(hzoom,'ActionPostCallback',@move_pos);

% apperance adjustment 
box on;
set(gcf, 'InvertHardCopy', 'off');

% appear
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
set_x_text(v0);
end

% ---------------------------------------------------------------------------
function modify_amplitude(hObject, eventdata)
v0= get(hObject,'Value');
v1=get(hObject,'min'); v2=get(hObject,'max');
Ly=(v2-v1)/2/5;

n=get(lbh,'Value');

if ( get(rbh_tsk_a,'Value')==get(rbh_tsk_a,'Max') )
   if v0>=0, fct=1+v0/Ly; else, fct=1+v0; end
   if fct==0 
      set(pid(n),'visible','off');
   else
      set(pid(n),'visible','on');
      set(pid(n),'ydata',(get(pid(n),'ydata')-Vpos(n))/Vamp(n)*fct+Vpos(n));
      Vamp(n)=fct;
   end
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
set_x_text(Vtime(n));
end

% ---------------------------------------------------------------------------
function set_x_text(v0)
if v0>get(hslider_x,'Min') & v0<get(hslider_x,'Max')
   set(htext_x,'string',num2str(v0),'fontsize',8, ...
   'position',[0.075+v0*0.83*(1-0.1-0.01)/50 0.001 0.04 0.02]);
   set(htext_x,'visible','on');
else
   set(htext_x,'visible','off');
end
end

% ---------------------------------------------------------------------------
function move_pos(obj, evd)
move_phase(obj,evd);
move_legend(obj,evd);
end

% ---------------------------------------------------------------------------
function [hlegd_line,hlegd_text]=plot_legend

set(gca,'ylim',get(gca,'ylim')*1.05);
scl_ylim=get(gca,'ylim'); scl_xlim=get(gca,'xlim');
Lx=(scl_xlim(2)-scl_xlim(1))/100; Ly=(scl_ylim(2)-scl_ylim(1))/11;

y0=scl_ylim(2)-Ly;
x1=2*Lx; x2=x1+2*Lx; x0=x2+Lx;

for n=1:nVar
    hlegd_line(n)=plot([x1 x2],[y0 y0]-(n-1)*Ly,V{n}.spec_line);
    hlegd_text(n)=text(x0,y0-(n-1)*Ly,[V{n}.name '(' num2str(V{n}.v0,'%3.2e') ')']);
end

end

% ---------------------------------------------------------------------------
function move_legend(obj, evd)

scl_ylim=get(gca,'ylim'); scl_xlim=get(gca,'xlim');
Lx=(scl_xlim(2)-scl_xlim(1))/100; Ly=(scl_ylim(2)-scl_ylim(1))/11;

y0=scl_ylim(2)-Ly;
x1=scl_xlim(1)+2*Lx; x2=x1+2*Lx; x0=x2+Lx;

for n=1:nVar
    set(hlegd_line(n),'xdata',[x1 x2],'ydata',[y0 y0]-(n-1)*Ly);
    set(hlegd_text(n),'position',[x0,y0-(n-1)*Ly,0]);
end

end

% ---------------------------------------------------------------------------
function [phspid,phstid]=plot_phase

if ~flag_phase, return; end

% borrowed from p1.m in saclab
  ylimm=get(gca,'YLim'); xlimm=get(gca,'XLim');
  ymin=ylimm(1,1); xmin=xlimm(1,1);
  ymax=ylimm(1,2); xmax=xlimm(1,2);
  yoff=.02*(ymax-ymin); xoff=.005*(xmax-xmin); xoff=0;

%phs_color=[0,0.75,0.75];
c_spec=colormap('lines'); nspec=size(c_spec,1);
NPHASE=numel(P);
for n=1:NPHASE
    phs_color=c_spec(mod(n-1,nspec)+1,:);
    ot=P(n).t; phase=P(n).name; fct=mod(n,6);
    hold on;
    phspid(n)=plot([P(n).t P(n).t],[ymin ymax],':','Color',phs_color,'Linewidth',1);
    phstid(n)=text(ot+xoff, ymin+yoff*fct, phase,'Color',phs_color,'FontName','times', ...
       'VerticalAlignment','bottom','FontSize',8);
    set(phspid(n),'ButtonDownFcn',{@remove_phase,phstid(n)});
end

end

% ---------------------------------------------------------------------------
function remove_phase(hObject, eventdata,ht)
set(hObject,'visible','off');
set(ht,'visible','off');
end

% ---------------------------------------------------------------------------
function move_phase(obj, evd)

if ~flag_phase, return; end

  ylimm=get(gca,'YLim'); xlimm=get(gca,'XLim');
  ymin=ylimm(1,1); xmin=xlimm(1,1);
  ymax=ylimm(1,2); xmax=xlimm(1,2);
  yoff=.02*(ymax-ymin); xoff=0;

NPHASE=numel(P);
for n=1:NPHASE
    fct=mod(n,6);
    set(phspid(n),'ydata',ylimm);
    pos=get(phstid(n),'position'); pos(2)=ymin+yoff*fct;
    set(phstid(n),'position',pos);
end

end

% ---------------------------------------------------------------------------

end
