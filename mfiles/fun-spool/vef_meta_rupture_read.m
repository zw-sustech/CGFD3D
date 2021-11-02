function [rxs,rys]=vef_meta_rupture_read(fault_fnm)

rx=nc_varget(fault_fnm,'x');
ry=nc_varget(fault_fnm,'y');
%rz=nc_varget(fault_fnm,'z');

rxs=[rx(1,:)';rx(:,end);rx(end,end:-1:1)';rx(end:-1:1,1)];
rys=[ry(1,:)';ry(:,end);ry(end,end:-1:1)';ry(end:-1:1,1)];

