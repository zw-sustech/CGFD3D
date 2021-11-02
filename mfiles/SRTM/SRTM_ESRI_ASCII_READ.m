function [ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value, topo] = SRTM_ESRI_ASCII_READ(filename) 
% SrtmEsriASCII-Process: A set of MATLAB program for processing Shuttle Radar
%   Topographic Mission (SRTM) topography data in Esri ASCII format.
%
%   (C) 2018 Wenzhong Cao, caowz@mail.ustc.edu.cn
%


% SRTM data can be download from http://srtm.csi.cgiar.org/.
%
%   SRTM data have two format:
%       1. GeoTiff
%       2. Esri ASCII
%   This program is used to do process for Esri ASCII format data. 
%
%   This function is use to read Esri ASCII format data:
%       ncols       : Number of columns
%       nrows       : Number of rows
%       xllcorner   : Longitude corresponding to the left point of the last row.
%       yllcorner   : Latitude  corresponding to the left point of the last row.
%       cellsize    : The width of a cell.
%       NODATA_value: This value indicates that the corresponding point has no data.
%       topo        : Topography data (2D array).
%   
    

HERE = 'SRTM_ESRI_ASCII_READ';

%% init variable.
ncols           = -1;
nrows           = -1;
xllcorner       = -1;
yllcorner       = -1;
cellsize        = -1;
NODATA_value    = -1;
topo            = [];

%% open file and read data.
fid = fopen(filename, 'r');
if fid<0
    disp(['Cannot open ' filename ' in function ' HERE]);
    return;
end

% read ncols/nrows.
[IntFormat]     = textscan(fid, '%s %d', 2);

if strcmp(cell2mat(IntFormat{1}(1)), 'ncols')
    ncols = IntFormat{2}(1);
else
    disp(['IntFormat error in function ' HERE '!!!']);
    return;
end

if strcmp(cell2mat(IntFormat{1}(2)), 'nrows')
    nrows = IntFormat{2}(2);
else
    disp(['IntFormat error in function ' HERE '!!!']);
    return;
end

% read xllcorner/yllcorner/cellsize/NODATA_value.
[FloatFormat] = textscan(fid, '%s %f', 3);

if strcmp(cell2mat(FloatFormat{1}(1)), 'xllcorner')
    xllcorner = FloatFormat{2}(1);
else
    disp(['FloatFormat error in function ' HERE '!!!']);
    return;
end

if strcmp(cell2mat(FloatFormat{1}(2)), 'yllcorner')
    yllcorner = FloatFormat{2}(2);
else
    disp(['FloatFormat error in function ' HERE '!!!']);
    return;
end

if strcmp(cell2mat(FloatFormat{1}(3)), 'cellsize')
    cellsize = FloatFormat{2}(3);
else
    disp(['FloatFormat error in function ' HERE '!!!']);
    return;
end

[IntFormat] = textscan(fid, '%s %d', 1);

if strcmp(cell2mat(IntFormat{1}(1)), 'NODATA_value')
    NODATA_value = IntFormat{2}(1);
else
    disp(['FloatFormat error in function ' HERE '!!!']);
    return;
end

readmat = repmat('%d', 1, ncols);
datamat = textscan(fid, readmat, nrows);

topo = cell2mat(datamat);

fclose(fid);

disp(['FINISHED @ ' HERE]);

