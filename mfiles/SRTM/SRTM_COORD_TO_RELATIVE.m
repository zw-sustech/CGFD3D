function [xb, yb, xe, ye, X1, Y1] = SRTM_COORD_TO_RELATIVE(x0, y0, X, Y, flag, scale)

    HERE='SRTM_COORD_TO_RELATIVE';

    X1 = X-x0;
    Y1 = Y-y0;

    xb = X1(1,1);
    xe = X1(1,end);

    yb = Y1(end,1);
    ye = Y1(1,1);

    disp(['FINISHED @ ' HERE]);

end
