function [cid,varargout]=vef_colorbar(hid,varargin)

args=varargin; nargs=nargin-1;

if isempty(hid), cid=colorbar; else cid=hid; end

% unpack input argurments
n=1;
while n<=nargs
  if ~ischar(args{n}), n=n+1; continue; end
  switch args{n}
  case {'location','loc'}
    set(cid,'location',args{n+1}); n=n+1;
  case {'position','pos'}
    set(cid,'position',args{n+1}); n=n+1;
  case {'xlim','ylim'}
    set(cid,args{n},args{n+1}); n=n+1;
  case {'xtick','ytick'}
    ytkloc=args{n+1};
    if strcmp(get(cid,[args{n}(1) 'limmode']),'manual')
       v0=get(cid,[args{n}(1) 'lim']); v0=v0(end);
       ytkloc(find(ytkloc>v0))=[];
       if v0-ytkloc(end)<(ytkloc(end)-ytkloc(end-1))/2, ytkloc(end)=[]; end
       ytkloc=[ytkloc v0]; ytkloc2=ytkloc(1:2:end); ytkloc2(end)=v0;
       set(cid,args{n},ytkloc);
       ytklab=get(cid,[args{n} 'label']);
    else
       set(cid,args{n},ytkloc);
    end
  case 'unit'
    id0=get(gcf,'CurrentAxes');
    set(gcf,'CurrentAxes',cid);
    pos=get(cid,'position');
    if pos(4)>pos(3)
       y0=get(cid,'ylim'); y0=(y0(end)-y0(1))/20+y0(end);
       x0=get(cid,'xlim'); x0=x0(end);
       tid=text(x0,y0,args{n+1});
    else
       y0=get(cid,'ylim'); y0=(y0(end)-y0(1))/2;
       x0=get(cid,'xlim'); x0=(x0(end)-x0(1))/20+x0(end);
       tid=text(x0,y0,args{n+1});
    end
    set(gcf,'CurrentAxes',id0);
    %ytlabel=get(cid,'yticklabel'); ytlabel(end,:)=[ytlabel(end,:) '(km/s^2)'];
  otherwise
  end
  n=n+1;
end

%keyboard;

nout=nargout-1;
if nout>=2
   varargout{1}=ytkloc; varargout{2}=ytklab;
end
if nout>=3
   varargout{3}=ytkloc2;
end
