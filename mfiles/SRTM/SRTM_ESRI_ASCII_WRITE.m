function SRTM_ESRI_ASCII_WRITE(fnm, ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value, topo)

    HERE='SRTM_ESRI_ASCII_WRITE';

    fid = fopen(fnm, 'w');
    if(fid<0)
        error(['Cannot open ' fnm ' @ ' HERE]);
    end

    fprintf(fid, 'ncols         %d\n', ncols);
    fprintf(fid, 'nrows         %d\n', nrows);
    fprintf(fid, 'xllcorner     %.12f\n', xllcorner);
    fprintf(fid, 'yllcorner     %.12f\n', yllcorner);
    fprintf(fid, 'cellsize      %.17f\n', cellsize);
    fprintf(fid, 'NODATA_value  %d\n', NODATA_value);

    for j=1:nrows
        for i=1:ncols
            fprintf(fid, '%d ', topo(j,i));
        end
        fprintf(fid, '\n');
    end

    fclose(fid);

    disp(['FINISHED @ ' HERE]);

end
