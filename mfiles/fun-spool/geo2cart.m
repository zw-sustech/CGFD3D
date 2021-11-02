function [x,y]=geo2cart(lat,lon,lat0,lon0,alpha,x0,y0)
% convert to Cartesian
% input unit degree
% output unit m

a=alpha/180*pi;

D=111319.5;

x=D*(cos(a)*(lat-lat0)+sin(a)*(lon-lon0).*cos(lat/180*pi));
y=D*(sin(a)*(lat-lat0)-cos(a)*(lon-lon0).*cos(lat/180*pi));

x=x+x0;
y=y+y0;
