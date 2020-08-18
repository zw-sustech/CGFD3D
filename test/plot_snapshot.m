
output_dir='/home/zhangw/work/cgfd_init_run/output';

snap_no = 2;
varname = 'Vz';

for it = 150 %0:10:500
  snap_fname = [output_dir, '/', ...
      'blk_1_snapshot_', num2str(snap_no), ...
      '_', varname, ...
      '_it', num2str(it), '.bin'];
  fid = fopen(snap_fname,'r');

  % read header
  % siz = fread(fid, 3, 'unit64');

  % read value
  var = fread(fid, inf, 'single');

  Vz = reshape(var, [100,100]);
  fread(fid);

  %surf(Vz);
  pcolor(Vz);
  title(['it=',num2str(it)]);
  %caxis([-1,1]*1e-25);
  colorbar
  drawnow;
end
