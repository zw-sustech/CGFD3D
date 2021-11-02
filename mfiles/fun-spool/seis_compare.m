function [varargout]=seis_compare(id,R,varargin);

%  Initial variable and flag
%S=varargin;
nVar=nargin-2;

if numel(id)>1
   error('under developing ...')
end

switch get(id,'type')
case 'root'
     scrsz = get(0,'ScreenSize');
     %fid=figure('position',scrsz, 'renderer','zbuffer','menubar','figure', ...
     %    'Name','Measurement of Records and Synthetic','NumberTitle','off');
     %fid=figure('position',[100 100 800 600], 'renderer','zbuffer','menubar','figure');
     fid=figure('renderer','zbuffer','menubar','figure');
     if isfield(R{1},'figtitle')
        set(fid,'Name',R{1}.figtitle);
     elseif isfield(R{1},'KEVNM')
        set(fid,'Name',['Comparison of record and synthetic for ', R{1}.KEVNM, '.', R{1}.KSTNM]);
     else
        set(fid,'Name','Comparison of record and synthetic');
     end
case 'figure'
     fid=id;
%case 'axes'
end

ntrace=numel(R);
elem_width=0.8; elem_height=0.8/ntrace; elem_left=0.1; elem_botm=0.1;
elem_pos=[elem_left,elem_botm,elem_width,elem_height];
incr_pos=[0,elem_height,0,0];

for n=1:ntrace
    aid(n)=axes('position',elem_pos+incr_pos*(n-1),'FontSize',8);
    hold on;
    pidr(n)=plot(R{n}.t,R{n}.v/R{n}.v0,R{n}.spec_line);
    for m=1:nVar
        S{m}=varargin{m}{n};
        pids{m}(n)=plot(S{m}.t,S{m}.v/S{m}.v0,S{m}.spec_line);
    end
    if exist('twin'); xlim(twin);end
    %legend(R{n}.name,S{n}.name);
    set(aid(n),'ButtonDownFcn',{@seis_measure,R{n},S});
    scl_ylim=get(gca,'ylim'); y0=scl_ylim(1)+(scl_ylim(2)-scl_ylim(1))/10;
    text(R{n}.t(1),y0,R{n}.name);
end
% disble ytick
set(aid,'yticklabel',[])
set(aid(2:end),'xticklabel',[])
set(aid,'box','on')

% xlabel
idxlab=get(aid(1),'XLabel'); set(idxlab,'string','Time (s)')

if isfield(R{1},'title')
   title(R{1}.title);
end

if nargout>0, varargout{1}=pidr; end
if nargout>1, varargout{2}=pids; end

%indx2=find(Tx>=scl_xlim2(1) & Tx<=scl_xlim2(2));
%
%NAX=numel(V);
%sp_color=colormap(lines);
%
%if 
%
%for n=1:NAX
%    hold on
%    lid(n)=plot(aid(n),V(n).t,V(n).v,spec_line);
%    xlim(scl_xlim1);
%end

