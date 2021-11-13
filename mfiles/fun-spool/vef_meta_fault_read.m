function [flat,flon]=vef_meta_fault_read(fault_list)

[flon,flat]=textread(fault_list,'%f %f');

