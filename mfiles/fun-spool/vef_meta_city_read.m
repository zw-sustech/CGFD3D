function [clat,clon,cnm]=vef_meta_city_read(city_list)

[clat,clon,cnm,cnmzh]=textread(city_list,'%f %f %s %s');

