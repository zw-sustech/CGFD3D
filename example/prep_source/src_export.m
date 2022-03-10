function src_export(fnm_ou, src)

% This function exports source model to .src file.
%  fnm_ou: string, the output file name;
%  src: source model. 
%     eventnm: string, event name,
%     number_of_source: number of point sources,
%     stf_is_discrete: flag, 0 or 1, is stf given by name or discrete values,
%     stf_time_length: time length of stf if stf_is_discrete == 0,
%     stf_dt: dt of discrete stf,
%     stf_nt: nt of discrete stf, 1 if stf_is_discrete == 0;
%     cmp_fi_mij: 1 force only, 2 stress only, 3 both force and stress,
%     mechansim_type: 0 by mij,  1 by strike/dip/rake and mu D A,
%     loc_coord_type: 0 grid index, 1 real coordinate,
%     loc_3dim: 0 meaning same as to loc_coord_type, 1 relative to free surface,
%     x,y,z_coord: [number_of_source] array, the location of each source;
%     t_start: [number_of_source] array, start time of each source;
%     stf_name: [number_of_source,:] array, stf name if stf_is_discrete == 0
%     stf_coefs: [10,number_of_source] array, stf coefs if stf_is_discrete == 0
%     Fx,Fy,Fz: [stf_nt,number_of_source] array, force values;
%     Mxx etc: [stf_nt,number_of_source] array, moment tensor values;
%     strike etc: [stf_nt,number_of_source] array, angles and Mu D A;
%       

%-- create output file
fid=fopen(fnm_ou,'w');

%-- 1st: event name
fprintf(fid, '%s\n', src.evtnm); 

%-- 2nd: number of source
fprintf(fid,'%d\n',src.number_of_source);

%-- 3rd: stf flags
if src.stf_is_discrete == 0
  fprintf(fid,'%d %g\n', 0, src.stf_time_length);
else
  fprintf(fid,'%d %g %d\n', 1, src.stf_dt, src.stf_nt);
end

%-- 4th: cmp and mechanism
fprintf(fid,'%d %d\n',src.cmp_fi_mij, src.mechansim_type);

%-- 5th: location type
fprintf(fid,'%d %d\n',src.loc_coord_type, src.loc_3dim);

%-- locations of each source
for i = 1 : src.number_of_source
    fprintf(fid,'%g %g %g\n',src.x_coord(i),src.y_coord(i),src.z_coord(i));
end

%-- cmp and mech of each source
if src.stf_is_discrete == 0

  for i = 1 : src.number_of_source
    fprintf(fid,'%g %s %g %g %g %g %g %g %g %g %g %g\n', ...
             src.t_start(i), src.stf_name(i,:),src.stf_coefs(:,i));

    if src.cmp_fi_mij ~= 2 % not stress only
        fprintf(fid,'%g %g %g', src.Fx(1,i), src.Fy(1,i), src.Fz(1,i));
    end

    if src.cmp_fi_mij >=2 
      if src.mechansim_type == 0
        fprintf(fid,' %g %g %g %g %g %g', ...
                src.Mxx(1,i),src.Myy(1,i),src.Mzz(1,i), ...
                src.Myz(1,i),src.Mxz(1,i),src.Mxy(1,i));
      else
        fprintf(fid,' %g %g %g %g %g %g', ...
                src.strike(1,i),src.dip(1,i),src.rake(1,i), ...
                src.mu(1,i),src.D(1,i),src.A(1,i));
      end
    end
              
    % return
    fprintf(fid, '\n');
  end

%- by discrete values
else

  for i = 1 : src.number_of_source
    fprintf(fid,'%g\n', src.t_start(i));

    for it = 1 : src.stf_nt

      if src.cmp_fi_mij ~= 2 % not stress only
          fprintf(fid,'%g %g %g', src.Fx(it,i), src.Fy(it,i), src.Fz(it,i));
      end

      if src.cmp_fi_mij >=2 
        if src.mechansim_type == 0
          fprintf(fid,' %g %g %g %g %g %g', ...
                  src.Mxx(it,i),src.Myy(it,i),src.Mzz(it,i), ...
                  src.Myz(it,i),src.Mxz(it,i),src.Mxy(it,i));
        else
          fprintf(fid,' %g %g %g %g %g %g', ...
                  src.strike(it,i),src.dip(it,i),src.rake(it,i), ...
                  src.mu(it,i),src.D(it,i),src.A(it,i));
        end
      end % stress
                
      % return
      fprintf(fid, '\n');

    end % stf_nt
  end % loop source

end % if name or discrete

fclose(fid);

end
