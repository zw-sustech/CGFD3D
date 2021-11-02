function [X, Y] = SRTM_COORD_CREATE(x0, y0, dx, dy, nx, ny)

    HERE='SRTM_COORD_CREATE';

    x1 = x0+(double(nx)-1)*dx;
    y1 = y0+(double(ny)-1)*dy;

    x = linspace(x0, x1, nx);
    y = fliplr(linspace(y0, y1, ny));

    for j=1:ny
        X(j, 1:nx) = x;
    end

    for i=1:nx
        Y(1:ny, i) = y;
    end

    disp(['FINISHED @ ' HERE]);

end
