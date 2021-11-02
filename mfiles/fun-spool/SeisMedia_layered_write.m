function SeisMedia_layered_write(L)

nx=length(L.x);
ny=length(L.y);
nlayer=size(L.z,3);

fid=fopen(L.filename,'w');

fprintf(fid,'%s\n','#######################################################################');
fprintf(fid,'%s\n','#                        model for layered structure                  #');
fprintf(fid,'%s\n','#######################################################################');

fprintf(fid,'distance2meter = %f\n', L.distance2meter);
fprintf(fid,'velocity2m/s     = %f\n', L.velocity2m    );
fprintf(fid,'density2kg/m^3   = %f\n', L.density2kg3   );

fprintf(fid,'number_of_layer = %i\n', nlayer ) ;
fprintf(fid,'horizontal_sampling = %i %i\n',nx,ny ) ;

fprintf(fid,'layer_meaning = %s\n',L.layermeaning ) ;
if L.withQs==1
   fprintf(fid,'QsF0 = %f\n', L.QsF0);
   fprintf(fid,'QsINF = %f\n',L.QsINF);
end

fprintf(fid,'%s\n','#######################################################################');
fprintf(fid,'%s\n','#                           media parameters                          #');
fprintf(fid,'%s\n','#######################################################################');
fprintf(fid,'%s\n','#  Vp(top) layer Vp(bottom)       Vs              density            Qs');
fprintf(fid,'%s\n','<anchor_media>');
for n=1:nlayer
    fprintf(fid,'%12.5g %-12.5g %12.5g %-12.5g %12.5g %-12.5g', ...
           L.Vp(n,:),L.Vs(n,:),L.rho(n,:));
    if L.withQs==1
       fprintf(fid,'%12.5g %-12.5g', L.Qs(n,:));
    end
    fprintf(fid,'\n');
end

fmt_topo='%16.7g %16.7g';
for n=1:nlayer; fmt_topo=[fmt_topo ' %16.7g']; end
fmt_topo=[fmt_topo '\n'];

fprintf(fid,'%s\n','#######################################################################');
fprintf(fid,'%s\n','#                        layer thickness or depth                     #');
fprintf(fid,'%s\n','#######################################################################');
fprintf(fid,'%s\n','<anchor_layer>');
for j=1:ny
    for i=1:nx
        fprintf(fid,fmt_topo,L.x(i),L.y(j),L.z(i,j,:));
    end
end
fclose(fid);

