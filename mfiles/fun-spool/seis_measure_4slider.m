function [varargout]=seis_measure(src,event,R,S,varargin)
% SIMPLE_GUI Select a data set from the pop-up menu, then
% click one of the plot-type push buttons. Clicking the button
% plots the selected data in the axes

%  Initialize and hide the GUI as it is being constructed.
% f = figure('Visible','off','Position',[360,500,450,285]);
R.name=strtrim(R.name); S.name=strtrim(S.name);

fid=figure('position',[700 100 800 600], 'renderer','zbuffer', ...
     'toolbar','figure','menubar','figure');
%set(fid,'Name',['Measurement of Records ' R.name ' and Synthetic ' S.name]);
set(fid,'Name',[R.name ' and ' S.name ' for ' R.KEVNM]);

%set(fid,'visible','off')

%htext  = uicontrol('Style','text','String','Select Data',...
%           'Position',[325,90,60,15]);

%ha = axes('Units','pixels','Position',[50,50,700,500]);
%ha = axes('Units','pixels');
ha = axes;

if 1
pidr=plot(R.t,R.v/R.v0,R.spec_line);
hold on
pids=plot(S.t,S.v/S.v0,S.spec_line);

scl_ylim=get(gca,'ylim'); scl_xlim=get(gca,'xlim');
Lx=(scl_xlim(2)-scl_xlim(1))/100; Ly=(scl_ylim(2)-scl_ylim(1))/11;

set(gca,'ylim',[scl_ylim(1)-Ly,scl_ylim(2)+Ly]);

y0=scl_ylim(2); y1=scl_ylim(2)-Ly;
x1=R.t(1)+2*Lx; x2=x1+2*Lx; x0=x2+Lx;

plot([x1 x2],[y0 y0],R.spec_line);
text(x0,y0,[R.name '(' num2str(R.v0,'%3.2e') ')']);

plot([x1 x2],[y1 y1],S.spec_line);
text(x0,y1,[S.name '(' num2str(S.v0,'%3.2e') ')']);

if isfield(R,'title')
   title(R.title);
end

if 1
hslider_r_time = uicontrol(fid,'Style','slider',...
                      'Max',20,'Min',0,'Value',0, ...
                      'SliderStep',[0.05 0.2], ...
                      'Unit','normalized','Position',[0.13 0.02 0.37 0.03], ...
                      'Callback',{@shift_time_axis,pidr,R.t});
hslider_s_time = uicontrol(fid,'Style','slider',...
                      'Max',20,'Min',0,'Value',0, ...
                      'SliderStep',[0.05 0.2], ...
                      'Unit','normalized','Position',[0.53 0.02 0.37 0.03], ...
                      'Callback',{@shift_time_axis,pids,S.t});

hslider_r_amp = uicontrol(fid,'Style','slider',...
                      'Max',1,'Min',-1,'Value',0, ...
                      'SliderStep',[0.01 0.1], ...
                      'Unit','normalized','Position',[0.95 0.55 0.02 0.38], ...
                      'Callback',{@modify_amplitude,pidr,R.v/R.v0});
hslider_s_amp = uicontrol(fid,'Style','slider',...
                      'Max',1,'Min',-1,'Value',0, ...
                      'SliderStep',[0.01 0.1], ...
                      'Unit','normalized','Position',[0.95 0.11 0.02 0.38], ...
                      'Callback',{@modify_amplitude,pids,S.v/S.v0});
end

if nargout>0, varargout{1}=pidr; end
if nargout>1, varargout{2}=pids; end
end

%end

% ---------------------------------------------------------------------------

function shift_time_axis(hObject, eventdata, lid, t)
v0= get(hObject,'Value');
set(lid,'xdata',t+v0);

function modify_amplitude(hObject, eventdata, lid, v)
v0= get(hObject,'Value');
v1=get(hObject,'min'); v2=get(hObject,'max');
Ly=(v2-v1)/2/5;
%y=get(lid,'ydata');
if v0>=0
   set(lid,'ydata',v*(1+v0/Ly));
else
   set(lid,'ydata',v*(1+v0));
end
%set(lid,'ydata',y*(1+v0
%% Proceed with callback...
