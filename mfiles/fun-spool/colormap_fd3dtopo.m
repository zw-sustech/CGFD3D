function cMap=colormap_fd3dtopo(varargin)
%
% colormap_fd3dtopo: Generate colormap.
%
% Usage: cMap=colormap_fd3dtopo(colorname)
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

% white-gray
% light blue
% green
% purple
% yellow-orange-red

nstep=10;

cMap=[];

cMap(end+1,:) = colorname_ImageMagick('snow1');
%%cMap(end+1,:) = colorname_ImageMagick('snow2');
%cMap(end+1,:) = colorname_ImageMagick('LightGray');
%
%cMap(end+1,:) = colorname_ImageMagick('SlateGray2');
%
%%cMap(end+1,:) = colorname_ImageMagick('PaleTurquoise1');
%%cMap(end+1,:) = colorname_ImageMagick('LightBlue1');
%%cMap(end+1,:) = colorname_ImageMagick('LightCyan2');
%
%cMap(end+1,:) = colorname_ImageMagick('LightSteelBlue');
%cMap(end+1,:) = colorname_ImageMagick('LightSkyBlue3');
%
%%cMap(end+1,:) = colorname_ImageMagick('SlateBlue');
%%cMap(end+1,:) = colorname_ImageMagick('MediumSlateBlue');
%cMap(end+1,:) = colorname_ImageMagick('LightSlateBlue');
%
%cMap(end+1,:) = colorname_ImageMagick('SteelBlue1');
%
%cMap(end+1,:) = colorname_ImageMagick('MediumTurquoise');

cMap(end+1,:) = colorname_ImageMagick('MediumSpringGreen');
cMap(end+1,:) = colorname_ImageMagick('SeaGreen2');

cMap(end+1,:) = colorname_ImageMagick('SpringGreen2');
cMap(end+1,:) = colorname_ImageMagick('LimeGreen');

%cMap(end+1,:) = colorname_ImageMagick('chartreuse1');
cMap(end+1,:) = colorname_ImageMagick('chartreuse2');
%cMap(end+1,:) = colorname_ImageMagick('chartreuse3');
cMap(end+1,:) = colorname_ImageMagick('LawnGreen');
%cMap(end+1,:) = colorname_ImageMagick('YellowGreen');
cMap(end+1,:) = colorname_ImageMagick('GreenYellow');

%cMap(end+1,:) = colorname_ImageMagick('yellow3');
%cMap(end+1,:) = colorname_ImageMagick('yellow2');

cMap(end+1,:) = colorname_ImageMagick('yellow');
cMap(end+1,:) = colorname_ImageMagick('gold');
cMap(end+1,:) = colorname_ImageMagick('DarkGoldenrod1');
cMap(end+1,:) = colorname_ImageMagick('orange');
cMap(end+1,:) = colorname_ImageMagick('DarkOrange');
cMap(end+1,:) = colorname_ImageMagick('OrangeRed');
%cMap(end+1,:) = colorname_ImageMagick('tomato');

cMap(end+1,:) = colorname_ImageMagick('red');
%cMap(end+1,:) = colorname_ImageMagick('OrangeRed2');
cMap(end+1,:) = colorname_ImageMagick('red2');
cMap(end+1,:) = colorname_ImageMagick('red3');
cMap(end+1,:) = colorname_ImageMagick('DarkRed');
