clc
clear all;
[ncols1, nrows1, xll1, yll1, cellsize, NODATA_value, topo1] = ...
    SRTM_ESRI_ASCII_READ('1.asc');

[ncols2, nrows2, xll2, yll2, cellsize, NODATA_value, topo2] = ...
    SRTM_ESRI_ASCII_READ('2.asc');

[ncols3, nrows3, xll3, yll3, cellsize, NODATA_value, topo3] = ...
    SRTM_ESRI_ASCII_READ('3.asc');

[ncols4, nrows4, xll4, yll4, cellsize, NODATA_value, topo4] = ...
    SRTM_ESRI_ASCII_READ('4.asc');
%%

topo=[topo1 topo2;topo3 topo4];
nrows=nrows1+nrows3;
ncols=ncols1+ncols2;
xll=xll1;
yll=yll3;


[X, Y] = SRTM_COORD_CREATE(xll, yll, cellsize, cellsize, ncols, nrows);

[xb, yb, xe, ye, X, Y] = SRTM_COORD_TO_RELATIVE(103, 28, X, Y);

figure(1)
surf(X, Y, topo);
shading interp;
xlim([xb xe]);
ylim([yb ye]);


% write data.
SRTM_ESRI_ASCII_WRITE('Yibin.asc', nrows, ncols, xll, yll, cellsize, NODATA_value, topo);