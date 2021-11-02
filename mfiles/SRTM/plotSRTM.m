clear all
clc

fid = 0;
%xll = 0;
%yll = 0;
%cellsize = 1;
%ncols = 6001;
%nrows = 6001;

% %% ----------------------------------------------------------------- %%
%% cut data%%
%% ----------------------------------------------------------------- %%
% read asc.
% [ncols, nrows, xll, yll, cellsize, NODATA_value, topo] = ...
% SRTM_ESRI_ASCII_READ('srtm_58_07.asc');
% 
% % cut data.
% [nx, ny, xb, yb, xe, ye, Topo] = ...
% SRTM_AREA_CUT(105, 28, 106, 30, cellsize, cellsize, xll, yll, ncols, nrows, topo);
% 
% % write data.
% SRTM_ESRI_ASCII_WRITE('4.asc', nx, ny, xb, yb, cellsize, NODATA_value, Topo);
% 
% % create coordinate.
% [X, Y] = SRTM_COORD_CREATE(xb, yb, cellsize, cellsize, nx, ny);
% 
% % trans coordinate to relative coordinate.
% [xb, yb, xe, ye, X, Y] = SRTM_COORD_TO_RELATIVE(103.5, 33, X, Y);
% 
% figure(1)
% surf(X, Y, Topo);
% axis equal;
% shading interp;
% xlim([xb xe]);
% %ylim([yb ye]);


%% ----------------------------------------------------------------- %%
% read & plot
% ----------------------------------------------------------------- %%
[ncols, nrows, xll, yll, cellsize, NODATA_value, topo] = ...
    SRTM_ESRI_ASCII_READ('Yibin.asc');

[X, Y] = SRTM_COORD_CREATE(xll, yll, cellsize, cellsize, ncols, nrows);

[xb, yb, xe, ye, X, Y] = SRTM_COORD_TO_RELATIVE(103.5, 33, X, Y);

fid = fid+1;
figure(fid)
surf(X, Y, topo);
shading interp;
xlim([xb xe]);
ylim([yb ye]);


