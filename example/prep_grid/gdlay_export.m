function gdlay_export(gdlay_file, ...
                      num_of_interfaces, ...
                      nx, ...
                      ny, ...
                      num_of_cell_per_layer, ...
                      dz_is_equal_of_layer, ...
                      x3d, ...
                      y3d, ...
                      z3d)

% This function exports gdlay to file

fid=fopen(gdlay_file,'w'); % Output file name 

%-- first line: how many interfaces
fprintf(fid,'%6d\n',num_of_interfaces);

%-- second line: how many cells
for ilay = 1 : num_of_interfaces - 1
  fprintf(fid,' %6d',num_of_cell_per_layer(ilay) );
end
fprintf(fid,'\n');

%-- third line: is dz equal of each layer
for ilay = 1 : num_of_interfaces - 1
  fprintf(fid,' %6d',dz_is_equal_of_layer(ilay) );
end
fprintf(fid,'\n');

%-- 4th line: nx, ny
fprintf(fid,'%6d %6d\n',nx, ny);

%-- others: z of interfaces
for n = 1 : num_of_interfaces
for j = 1 : ny
for i = 1 : nx
    %-- seems the first layer is the deepest one
    fprintf(fid,'%12.2f %12.2f %12.2f\n', ...
        x3d(i,j,n), y3d(i,j,n),z3d(i,j,n));
end
end
end

fclose(fid);

end % function
