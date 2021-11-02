function [cMap,ticks,nelems]=colormap_stepseq(varargin)
%
% colormap_stepseq: Generate colormap.
%
% Usage: cMap=colormap_stepseq(colorname)
%

% Major ChangeLog:
%   2009-01-09 Wei Zhang
%     * Added help information, but uncomplete.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Date$
% $Revision$
% $LastChangedBy$
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default parameters
nstep=5; nseq=5;
nratio=[1,1,1];

% unpack input arguement
args=varargin; nargs=nargin; n=1;
while n<=nargs
  if ~ischar(args{n}), n=n+1; continue; end
  switch lower(args{n})
  case 'nstep'
    nstep=args{n+1}; n=n+1;
  case 'nseq'
    nseq=args{n+1}; n=n+1;
  case 'nratio'
    nratio=args{n+1}; n=n+1;
  end
  n=n+1;
end

nlen=length(nratio);
if nlen==3
   nratio(nstep)=nratio(nlen);
   nratio(2:nstep-1)=nratio(2);
else
   nratio(nlen+1:nstep)=nratio(nlen);
end
ticks=[0,cumsum(nratio)/min(min(nratio))];
nelems=sum(nratio)/min(min(nratio));
nratio=nratio/min(min(nratio))*nseq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        generate cMap                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cMap=[];
n=0;

% -- use http://www.wellstyled.com/tools/colorscheme2/index-en.html --

% purple
n=n+1; nseq=nratio(n);
c1=name2rgb('white'); c2=name2rgb('MediumPurple2');
%c1=name2rgb('white'); c2=hex2rgb('CC0099');
C=makeColorMap(c1,c2,nseq); cMap=[cMap; C];

% blue
if nstep>=7
n=n+1; nseq=nratio(n);
%c1=name2rgb('SkyBlue'); c2=name2rgb('RoyalBlue');
c1=hex2rgb('9191FF'); c2=hex2rgb('1919B3');
C=makeColorMap(c1,c2,nseq); cMap=[cMap; C];
end

% green 0
if nstep>=3
n=n+1; nseq=nratio(n);
c1=hex2rgb('B8FF4D'); c2=hex2rgb('6BB300');
C=makeColorMap(c1,c2,nseq); cMap=[cMap; C];
end

% green
if nstep>=9
n=n+1; nseq=nratio(n);
%c1=name2rgb('MediumSpringGreen'); c2=name2rgb('ForestGreen');
c1=hex2rgb('80FF80'); c2=hex2rgb('008F00');
C=makeColorMap(c1,c2,nseq); cMap=[cMap; C];
end

% yellow
if nstep>=5
n=n+1; nseq=nratio(n);
%c1=name2rgb('khaki2'); c2=name2rgb('yellow'); c3=name2rgb('yellow2');
%c1=hex2rgb('FFFF80'); c2=hex2rgb('FFFF00'); c3=hex2rgb('E6E600');
%c1=hex2rgb('FFFF73'); c2=hex2rgb('FFFF00'); c3=hex2rgb('E6E600');
%c1=hex2rgb('FFFF66'); c2=hex2rgb('FFFF00'); c3=hex2rgb('E6E600');
%c1=hex2rgb('FFFF66'); c2=hex2rgb('FFFF00'); c3=hex2rgb('CCCC00');
%C=makeColorMap(c1,c2,c3,nseq); cMap=[cMap; C];
c1=hex2rgb('FFFF73'); c2=hex2rgb('D9D900');
C=makeColorMap(c1,c2,nseq); cMap=[cMap; C];
end

% orange
if nstep>=4
n=n+1; nseq=nratio(n);
%c1=name2rgb('goldenrod1'); c2=name2rgb('DarkOrange');
%c1=hex2rgb('FFD980'); c2=hex2rgb('FFB200');
%c1=hex2rgb('FFDB4D'); c2=hex2rgb('D9AD00');
%c1=hex2rgb('FFDB4D'); c2=hex2rgb('BF9900');
%c1=hex2rgb('FFCD59'); c2=hex2rgb('B37D00');
%c1=hex2rgb('FFD166'); c2=hex2rgb('B37D00');
%c1=hex2rgb('FFA64D'); c2=hex2rgb('FF8000');
c1=hex2rgb('FFC266'); c2=hex2rgb('F27A00');
C=makeColorMap(c1,c2,nseq); cMap=[cMap; C];
%c1=hex2rgb('FFA64D'); c2=hex2rgb('FF8000'); c3=hex2rgb('B35A00');
%C=makeColorMap(c1,c2,c3,nseq); cMap=[cMap; C];
end

%% tomato
%if 0 %nstep>=9
%%c1=name2rgb('DarkOrange2'); c2=name2rgb('OrangeRed');
%%c1=hex2rgb('FFB380'); c2=hex2rgb('FF6600');
%c1=hex2rgb('FF944D'); c2=hex2rgb('D95700');
%C=makeColorMap(c1,c2,nseq); cMap=[cMap; C];
%end

% pink
if nstep>=10
n=n+1; nseq=nratio(n);
%c1=name2rgb('pink'); c2=name2rgb('DeepPink');
%c1=hex2rgb('FF59A3'); c2=hex2rgb('E60066');
c1=hex2rgb('FF80DF'); c2=hex2rgb('E60066');
C=makeColorMap(c1,c2,nseq); cMap=[cMap; C];
end

% red 0
if nstep>=11
n=n+1; nseq=nratio(n);
%c1=hex2rgb('FF704D'); c2=hex2rgb('FF3300');
%c1=hex2rgb('FF703D'); c2=hex2rgb('FF3d0D');
%c1=hex2rgb('FF8566'); c2=hex2rgb('FF3d0D');
c1=hex2rgb('FF9980'); c2=hex2rgb('F24E24');
C=makeColorMap(c1,c2,nseq); cMap=[cMap; C];
end

% red
n=n+1; nseq=nratio(n);
%c1=name2rgb('red'); c2=name2rgb('red2');
%c1=hex2rgb('FF4040'); c2=hex2rgb('FF0000');
c1=hex2rgb('FF5959'); c2=hex2rgb('FF0000');
C=makeColorMap(c1,c2,nseq); cMap=[cMap; C];

% dark red
if nstep>=6
n=n+1; nseq=nratio(n);
%c1=name2rgb('red3'); c2=name2rgb('DarkRed');
%c1=hex2rgb('F20000'); c2=hex2rgb('D90000');
%c1=hex2rgb('D94C4C'); c2=hex2rgb('B30000');
c1=hex2rgb('E66767'); c2=hex2rgb('D90000');
C=makeColorMap(c1,c2,nseq); cMap=[cMap; C];
end

% over dark red
if nstep>=8
n=n+1; nseq=nratio(n);
%c1=hex2rgb('8C2A2A'); c2=hex2rgb('730000');
c1=hex2rgb('B35050'); c2=hex2rgb('730000');
C=makeColorMap(c1,c2,nseq); cMap=[cMap; C];
end


