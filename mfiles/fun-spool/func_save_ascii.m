function [ierr] = func_save_ascii(T,V,fnm_ou)

ierr=0;

nt = length(T);

ncmp = size(V,2);

if nt ~= size(V,1)
    display('Error: dimension mismatch');
    length(T)
    size(V,1)
    size(V,2)
    ierr=-1;
    return;
end

if exist(fnm_ou,'file')
    display(['Warning: ', fnm_ou, ' exists, overwriting now ...']);
end

fid = fopen(fnm_ou,'w');

fmt_str = '%g';
for n=1:ncmp
    fmt_str = [ fmt_str ' %g'];
end
fmt_str = [ fmt_str '\n' ];

for n=1:nt
    fprintf(fid,fmt_str, T(n), V(n,1:ncmp));
end

fclose(fid);

disp(['finish func_save_ascii ', fnm_ou]);
