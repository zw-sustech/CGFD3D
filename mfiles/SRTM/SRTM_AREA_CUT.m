function [nx, ny, xb, yb, xe, ye, Topo] = SRTM_AREA_CUT(x0, y0, x1, y1, dx, dy, xll, yll, ncols, nrows, Ltopo)

    HERE='SRTM_AREA_CUT';

    if(x1<x0)
        disp(['ERROR:   Latitude right side is smaller than left side.']);
        return;
    end

    if(y1<y0)
        disp(['ERROR:   Longitude top side is smaller than bottom side.']);
        return;
    end

    xrr = xll+(double(ncols)-1)*dx;
    yrr = yll+(double(nrows)-1)*dy;


    if(x0<xll)
        x0=xll;
        disp(['WARRING: Latitude is out of range (left side.)!']);
    end

    if(x1>xrr)
        x1 = xrr;
        disp(['WARRING: Latitude is out of range (right side.)!']);
    end

    if(y0<yll)
        y0=yll;
        disp(['WARRING: Longitude is out of range (bottom side.)!']);
    end

    if(y1>yrr)
        y1 = yrr;
        disp(['WARRING: Longitude is out of range (top side.)!']);
    end

    ib = floor((x0-xll)/dx)+1;
    xb = xll + (double(ib)-1)*dx;

    jb = floor((y0-yll)/dy)+1;
    yb = yll + (double(jb)-1)*dy;
    
    ie = ib + floor((x1-x0)/dx);
    xe = xll + (double(ie)-1)*dx;

    je = jb + floor((y1-y0)/dy);
    ye = yll + (double(je)-1)*dy;

    nx = ie - ib + 1;
    ny = je - jb + 1;


    Topo = Ltopo(nrows-je+1:nrows-jb+1, ib:ie);

    disp(['FINISHED @ ' HERE]);

end
