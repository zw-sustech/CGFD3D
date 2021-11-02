function h=vef_figure(varargin)

% ----------------------- parameter -----------------------
args=varargin; nargs=numel(args);
n=1;
while n<=nargs

if ~ischar(args{n}), n=n+1; continue; end

switch args{n}
case 'bgcolor', bgcolor=args{n+1};
end

n=n+1;
end

% ----------------------- figure -----------------------
h=figure;
set(h,'BackingStore','on');
set(h,'renderer','zbuffer');
set(h,'menubar','none');
set(h,'toolbar','figure');
%set(h, 'PaperPositionMode', 'manual');
%set(h,'PaperUnits','points')
%set(h,'PaperPosition',[0 0 1024 768])

%set(h, 'PaperPositionMode', 'manual');
%set(h,'PaperUnits','inches')
%set(h,'PaperPosition',[0 0 3.3 2.5])

if exist('bgcolor','var')
   set(h,'color',bgcolor);
end
