function [fnm_out]=get_fnm_metric(pnm,n_i,n_j,n_k,varargin)

fnm_out=[pnm '/' 'coord_mpi'  ...
         num2str(n_i,'%2.2i') ...
         num2str(n_j,'%2.2i') ...
         num2str(n_k,'%2.2i') ...
         '.nc'];

nvin=nargin-4;
if nvin==0, return; end

% return old filename if 'R1' argue exits
n=1;
while n<=nvin

  if numel(varargin{n})~=1, n=n+1; continue; end
  if isnumeric(varargin{n}), n=n+1; continue; end
  
  switch varargin{n}
  case {'r1','R1'}
    fnm_out=[pnm '/' 'swmpi'  ...
             num2str(n_i,'%2.2i') ...
             num2str(n_j,'%2.2i') ...
             num2str(n_k,'%2.2i') ...
             '_coord.nc'];
  end

end
