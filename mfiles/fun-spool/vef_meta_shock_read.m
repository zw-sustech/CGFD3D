function [lat,lon,nm]=vef_meta_shock_read(shock_list)

[lon,lat,nm]=textread(shock_list,'%f %f %s');

