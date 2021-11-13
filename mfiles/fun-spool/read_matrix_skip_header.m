function A=read_matrix_skip_header(fnm,nskip,fmt_str,rsize);
% read ascii matrix with nskip header

fid=fopen(fnm,'r');
for n=1:nskip, fgetl(fid); end

A=fscanf(fid,fmt_str,rsize);

fclose(fid);
