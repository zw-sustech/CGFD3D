function img = imcapture( h, opt, dpi, opt2, opt3)
% IMCAPTURE do screen captures at controllable resolution using the undocumented "hardcopy" built-in function.
%
% USAGE:
%   IMG = IMCAPTURE(H) gets a screen capture from object H, where H is a handle to a figure, axis or an image
%                   When H is an axis handle and OPT ~= 'all' (see below) the capture is done on
%                   that axis alone, no matter how many more axis the Figure contains. 
%   IMG = IMCAPTURE or IMG = IMCAPTURE([]) operates on the current figure.
%   IMG = IMCAPTURE(H, OPT) selects one of three possible operations depending on the OPT string
%           OPT == 'img' returns an image corresponding to the axes contents and a size determined
%                   by the figure size.
%           OPT == 'imgAx' returns an image corresponding to the displyed image plus its labels.
%                   This option tries to fill a enormous gap in Matlab base which is the non availability of
%                   an easy way of geting a raster image that respects the data original aspect ratio.
%           OPT == 'all' returns an image with all the figure's contents and size determined by the figure size.
%   IMG = IMCAPTURE(H,OPT,DPI) do the screen capture at DPI resolution. DPI can be either a string or numeric
%                   If DPI is not provided, the screen capture is done at 150 dpi
%                   If DPI = 0, returns an image that has exactly the same width and height of the original CData
%                   You can use this, for example, to convert an indexed image to RGB.
%                   Note, this should be equivalent to a call to GETFRAME, but maybe it fails less in
%                   returning EXACTLY an image of the same size of CData (R13, as usual, is better than R14).
%   IMG = IMCAPTURE(H,'img',[MROWS NCOLS]) returns an image of the size specified by [mrows ncols]
%                   This a 10 times faster option to the use of imresize(IMG,[mrows ncols],method)
%                   when method is either 'bilinear' or 'bicubic'. Not to mention the memory consumption.
%
%   IMCAPTURE(H,OPT,FNAME) Save the screen capture on file FNAME using EPS, AI or EMF formats
%   IMCAPTURE(H,OPT,FNAME,DPI,SIZE) Save the screen capture on file FNAME using EPS, AI or EMF formats
%           OPT == 'img' or 'imgAx' their meaning is the same as explained above.
%           FNAME   is a file name with one of the following 4 possible extensions .ps .eps .ai .emf
%                   Obviously, the extension type selects the output format.
%           DPI     Same as above except that if not present, DPI defaults to 300.
%           SIZE    Is a two element vector [WIDTH HEIGHT] with the image dimensions in centimeters.
%                   If not provided, image dimensions are picked from the Figure's 'PaperPosition' property.
%
%   IMCAPTURE(H,OPT) Copy the figure as a Windows Enhanced Metafile and place it in the clipboard.
%                   This works only in Windows.
%
%   NOTE. This function also works with the 'imgAx' option for plots and surfaces but DPI settings
%   is more "fluid" because I cannot exactly determine what is the image.
%
%   EXAMPLE:
%       load clown
%       h = imagesc(X);
%       set(gcf,'colormap',map)
%       x = 27*cos(-pi:.1:pi) + 55;
%       y = 27*sin(-pi:.1:pi) + 120;
%       patch('XData',x,'YData',y,'FaceColor','y');     % Let's change the clown nose color
%       % Notice how the yellow circle as turned into an ellipse.
%       % Now do an image screen capture at 300 dpi
%       rgb=imcapture(h,'img',300);
%       % In order to confirm the correct aspect ratio you need either to call axis image
%       figure; image(rgb); axis image
%       % or save into an raster format an open it with outside Matlab. E.G
%       imwrite(rgb,'lixo.jpg')

%   HOW DOES IT WORKS?
%       A lot of testings revealed that the array size returned by 'hardcopy()' function
%       depends on the 'PaperPosition' property.
%       BTW the help hardcopy says
%           "HARDCOPY(H,'filename','format') saves the figure window with handle H
%            to the designated filename in the specified format."
%       and
%           "Do NOT use this function directly. Use PRINT instead."
%
%       Well, fortunately the first is false and I didn't follow the advise.
%       But let us return to the 'PaperPosition' property. By default figures are created
%       (at least in my system) with this values in centimeters
%       get(hFig,'PaperPosition') = [0.6345    6.3452   20.3046   15.2284]
%       Using the example above (the clown) X = 200x320
%           rgb= hardcopy(get(h,'parent'),'lixo.jpg','-dzbuffer');
%       shows that rgb = 900x1201x3. And how did hardcopy picked these dimensions?
%       We can find that X & Y resolution are different
%           DPIx = round(320 / 20.3046 * 2.54) = 40
%           DPIy = round(200 / 15.2284 * 2.54) = 33
%       Shit, we have different resolutions - what are these guys doing?
%       But see this
%           DPIx = round(1201 / 20.3046 * 2.54) = 150
%       So the capture is carried out at 150 dpi, and using this on the DPIy expression
%           N = 15.22842 * 150 / 2.54 = 899
%       We nearly got it. 899 instead of 900. This is may be due to inch<->centimeters
%       rounding error ... or to a bug. 
%       OK, now that we found how the size is picked, the trick is to play arround with
%       the 'PaperPosition' in order to get what we want.
%
%   AUTHOR
%       Joaquim Luis (jluis@ualg.pt)    12-Dez-2006
%       University of Algarve
%
%       Revisions
%           13-Dez-2006     Let it work also with plots
%           28-Dez-2006     Added vector graphics output and copy to Clipboard (Windows only)
%           12-Jan-2007     Improved "fixed size" option. Axes are set visible when their visibility was 'off'

    hFig = [];      hAxes = [];
    if (nargin == 0 || isempty(h)),     h = get(0,'CurrentFigure');    end
    if (~ishandle(h))
        error('imcapture:a','First argument is not a valid handle')
    end
    if (strcmp(get(h,'Type'),'figure'))
    	hAxes = findobj(h,'Type','axes');
    elseif (strcmp(get(h,'Type'),'axes'))
        hAxes = h;
        h = get(h,'Parent');
    elseif (strcmp(get(h,'Type'),'image'))
        hAxes = get(h,'Parent');
        h = get(hAxes,'Parent');
    else
        h = [];
    end
    if (~ishandle(h))
        error('imcapture:a','First argument is not a valid Fig/Axes/Image handle')
    end
    if (nargin <= 1),   opt = [];   end
    
    inputargs{1} = h;
    renderer = get(h, 'Renderer');
    
    if (nargout == 1)                   % Raster mode. We expect to return a RGB array
        inputargs{4} = '-r150';         % Default value for the case we got none in input
        inputargs{2} = 'lixo.jpg';      % The name doesn't really matter, but we need one.
        if strcmp(renderer,'painters')
            renderer = 'zbuffer';
        end
        inputargs{3} = ['-d' renderer];
        if (nargin == 3)
            if (isnumeric(dpi) && numel(dpi) == 1)
                inputargs{4} = ['-r' sprintf('%d',dpi)];
            elseif (isnumeric(dpi) && numel(dpi) == 2)      % New image size in [mrows ncols]
                inputargs{4} = dpi;
            elseif (ischar(dpi))
                inputargs{4} = ['-r' dpi];
            else
                error('imcapture:a','third argument must be a ONE or TWO elements vector, OR a char string')
            end
        end
    elseif (nargout == 0 && nargin == 2)    % Clipboard
        if (~strncmp(computer,'PC',2))
            error('imcapture:a','Copying to clipboard is only possible in windows')
        end
        inputargs{2} = '';
        inputargs{3} = '-dmeta';
        inputargs{4} = '-r300';     % I think that it realy doesn't matter in this case
    else                            % Vector graphics mode. Returns nothing but save capture on file
        if (nargin < 3),    error('imcapture:a','Missing output filename.');  end
        if (~ischar(dpi)),  error('imcapture:a','Third argument must be a string with the filename');  end
        inputargs{4} = '-r300';     % Default value for the case we got none in input
        [pato,fname,ext] = fileparts(dpi);
        switch lower(ext)
            case {'.ps','.eps'},    inputargs{3} = '-deps2';
            case '.ai',             inputargs{3} = '-dill';
            case '.emf',            inputargs{3} = '-dmeta';
            otherwise
                error('imcapture:a','Illegal extension in filename. Legal extensions are .ps, .eps, .ai or .emf')
        end
        inputargs{2} = dpi;
        if (nargin >= 4)
            if (isnumeric(opt2)),       inputargs{4} = ['-r' sprintf('%d',opt2)];
            else                        inputargs{4} = ['-r' opt2];
            end
        end
        if (nargin == 5)
            if (~isnumeric(opt3) || numel(opt3) ~= 2)
                error('imcapture:a','Fifth argument must be a 2 element vector with image width and height')
            end
            inputargs{5} = opt3;
        end
    end
    
    msg = [];
    if (numel(hAxes) == 1 && strcmp( get(hAxes,'Visible'),'off') && (nargin == 1 || isempty(opt)) )
        % Default to 'imgOnly' when we have only one image with invisible axes
        opt = 'img';
    end

    if (nargin > 1)
        switch opt
            case 'img'      % Capture the image only
                [img,msg] = imgOnly([],hAxes,inputargs{:});
            case 'imgAx'    % Capture an image including the Labels
                [img,msg] = imgOnly('nikles',hAxes,inputargs{:});
            case 'all'      % Capture everything in Figure
                img = allInFig(inputargs{:});
            otherwise
                msg = 'Second argument is not a recognized option';
        end
    else
        img = allInFig(inputargs{:});
    end
    
    if (~isempty(msg))      % If we had an error inside imgOnly()
        error('imcapture:a',msg);        img = [];
    end

% ------------------------------------------------------------------    
function [img, msg] = imgOnly(opt, hAxes, varargin)
    % Capture the image, and optionaly the frame, mantaining the original image aspect ratio.
    % We do that be messing with the Figure's 'PaperPosition' property
    h = varargin{1};    msg = [];
    if (isempty(hAxes) || numel(hAxes) > 1)
        msg = 'With the selected options the figure must contain one, and one ONLY axes';
        return
    end
    im = get(findobj(h,'Type','image'),'CData');
	if (~isempty(im))
        nx = size(im,2);                ny = size(im,1);
    else                    % We have something else. A plot, a surface, etc ...
        axUnit = get(hAxes,'Units');    set(hAxes,'Units','pixels')
        axPos = get(hAxes,'pos');       set(hAxes,'Units',axUnit)
        nx = axPos(3);                  ny = axPos(4);
	end

    PU = get(h,'paperunits');       set(h,'paperunits','inch')
    pp = get(h,'paperposition');    PPM = get(h,'PaperPositionMode');
    dpi = (nx / pp(3));             % I used to have a round() here but this pobably more correct
    % Here is the key point of all this manip.
    set(h,'paperposition',[pp(1:3) ny / dpi])

    axUnit = get(hAxes,'Units');
    axPos = get(hAxes,'pos');           % Save this because we will have to restore it later
    set(hAxes,'Units','Normalized')     % This is the default, but be sure
    fig_c = get(h,'Color');       set(h,'Color','w')
    
    all_axes = findobj(h,'Type','axes');
    % If there are more than one axis, put them out of sight so that they wont interfere
    % (Just make them invisible is not enough)
    if (numel(all_axes) > 1)
        other_axes = setxor(hAxes,all_axes);
        otherPos = get(other_axes,'Pos');
        bakAxPos = otherPos;
        if (~iscell(otherPos))
            otherPos(1) = 5;
            set(other_axes,'Pos',otherPos)
        else
            for (i=1:numel(other_axes))
                otherPos{i}(1) = 5;
                set(other_axes(i),'Pos',otherPos{i})
            end
        end
    end
    
    have_frames = false;
    axVis = get(hAxes,'Visible');
    if (isempty(opt))                   % Pure Image only capture. Even if axes are visible, ignore them
        set(hAxes,'pos',[0 0 1 1],'Visible','off')
    elseif (strcmp(get(hAxes,'Visible'),'on'))  % Try to capture an image that respects the data aspect ratio
        have_frames = true;
		h_Xlabel = get(hAxes,'Xlabel');         h_Ylabel = get(hAxes,'Ylabel');
		units_save = get(h_Xlabel,'units');
		set(h_Xlabel,'units','pixels');         set(h_Ylabel,'units','pixels');
		Xlabel_pos = get(h_Xlabel,'pos');
		Ylabel_pos = get(h_Ylabel,'Extent');
		
		if (abs(Ylabel_pos(1)) < 30)    % Stupid hack, but there is a bug somewhere
            Ylabel_pos(1) = 30;
		end
		
		y_margin = abs(Xlabel_pos(2))+get(h_Xlabel,'Margin');  % To hold the Xlabel height
		x_margin = abs(Ylabel_pos(1))+get(h_Ylabel,'Margin');  % To hold the Ylabel width
		y_margin = min(max(y_margin,20),30);            % Another hack due to the LabelPos non-sense
        
        figUnit = get(h,'Units');        set(h,'Units','pixels')
        figPos = get(h,'pos');           set(h,'Units',figUnit)
        x0 = x_margin / figPos(3);
        y0 = y_margin / figPos(4);
        set(hAxes,'pos',[x0 y0 1-[x0 y0]-1e-2])
        set(h_Xlabel,'units',units_save);     set(h_Ylabel,'units',units_save);
    else            % Dumb choice. 'imgAx' selected but axes are invisible. Default to Image only
        set(hAxes,'pos',[0 0 1 1],'Visible','off')
    end
    
    confirm = false;
    try
        if (strcmp(varargin{4},'-r0'))              % One-to-one capture
            varargin{4} = ['-r' sprintf('%d',round(dpi))];
            confirm = true;
            mrows = ny;            ncols = nx;      % To use in "confirm"
        elseif (numel(varargin{4}) == 2)            % New size in mrows ncols
            mrows = varargin{4}(1);
            ncols = varargin{4}(2);
            if (~have_frames)
                set(h,'paperposition',[pp(1:2) ncols/dpi mrows/dpi])
                confirm = true;
            else                        % This is kind of very idiot selection, but let it go
                wdt = pp(3);          hgt = pp(3) * mrows/ncols;
                set(h,'paperposition',[pp(1:2) wdt hgt])
            end
            varargin{4} = ['-r' sprintf('%d',round(dpi))];
        else                            % Fourth arg contains the dpi
            if (have_frames)
                wdt = pp(3);          hgt = pp(3) * ny/nx;
                set(h,'paperposition',[pp(1:2) wdt hgt])
            end
        end
        if (numel(varargin) == 5)       % Vector graphics formats
            set(h,'paperposition',[pp(1:2) varargin{5}*2.54])       % I warned in doc to use CM dimensions
            varargin(5) = [];           % We don't want it to go into hardcopy
        end
        
        img = hardcopy( varargin{:} );      % Capture
        if (confirm)                        % We asked for a pre-determined size. Check that the result is correct
            dy = mrows - size(img,1);       % DX & DY should be zero or one (when it buggs).
            dx = ncols - size(img,2);
            if (dx ~= 0 || dy ~= 0)         % ML failed (probably R14). Repeat to make it obey
                mrows_desBUG = mrows + dy;
                ncols_desBUG = ncols + dx;
                set(h,'paperposition',[pp(1:2) ncols_desBUG/dpi mrows_desBUG/dpi])
                img = hardcopy( varargin{:} );      % Insist
            end
        end
    catch                                   % If it screws, restore original Fig properties anyway
        set(hAxes,'Units',axUnit,'pos',axPos,'Visible','on')
        set(h,'paperposition',pp,'paperunits',PU,'PaperPositionMode',PPM,'Color',fig_c)
        msg = lasterr;      img = [];
    end
   
    % If there are more than one axis, bring them to their's original positions
    if (numel(all_axes) > 1)
        if (~iscell(otherPos))
            set(other_axes,'Pos',bakAxPos)
        else
            for (i=1:numel(other_axes)),    set(other_axes(i),'Pos',bakAxPos{i});  end
        end
    end
    
    % Reset the original fig properties
    set(hAxes,'Units',axUnit,'pos',axPos,'Visible',axVis)
    set(h,'paperposition',pp,'paperunits',PU,'PaperPositionMode',PPM,'Color',fig_c)
    
% ------------------------------------------------------------------    
function img = allInFig(varargin)
    % Get everything in the Figure
    h = varargin{1};
    fig_c = get(h,'Color');       set(h,'Color','w')
    if (numel(varargin) == 3)
        varargin{4} = '-r150';
    end
    img = hardcopy( varargin{:} );    
    set(h,'Color',fig_c)
