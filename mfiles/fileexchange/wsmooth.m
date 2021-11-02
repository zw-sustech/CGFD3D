function [zs,L] = wsmooth(z,x,y,L)

%WSMOOTH   1D and 2D robust smoothing.
%   ZS = WSMOOTH(Z,X,L) smoothes the signal Z(X) using the smoothing
%   parameter L. Z and X must be vector arrays of same length. L must be a
%   negative or positive real scalar (the program actually uses 10^L). The
%   larger L is, the smoother the output will be. Typical values for L are
%   in the range [0,10]. 
%
%   ZS = WSMOOTH(Z,X,Y,L) smoothes the surface Z(X,Y). X and Y must be
%   vector arrays whose lengths match the size of matrix Z, i.e.:
%   LENGTH(X) = SIZE(Z,2) and LENGTH(Y) = SIZE(Z,1)
%
%   ZS = WSMOOTH(Z,L) uses unit spacings for X and/or Y.
%
%   ZS = WSMOOTH(Z,X,Y) or ZS = WSMOOTH(Z,X) or ZS = WSMOOTH(Z): if L is
%   omitted, it is automatically determined using the generalized cross-
%   validation (GCV) method (see details in the program, line 260). Because
%   eigenvalues are calculated, this may take a few seconds with large
%   arrays. [ZS,L] = WSMOOTH(...) also returns the calculated value for L
%   so that you can fine-tune the smoothing subsequently if wanted.
%   
%   Notes:
%   -----
%   WSMOOTH also works if some data are missing (see examples); data are
%   considered missing if they are not finite (NaN or Inf). WSMOOTH can
%   thus be used as an alternative to INPAINT_NANS (a John D'Errico's
%   function): Z = INPAINT_NANS(ZNAN) and Z = WSMOOTH(ZNAN,L) with small L
%   give nearly similar results. Smoothing level, however, can be adjusted
%   with WSMOOTH if necessary.
%
%   This program is based on the Whittaker's smoother. The 1D algorithm has
%   been well described by Eilers in Anal. Chem. 2003, 75, 3631-3636. This
%   algorithm has been here adapted to 2D. Note that a 3D update could be
%   easily written (but beware the "out of memory" errors).
%   
%   Examples:
%   --------
%   % 1-D example
%   x = 0:100;
%   y = cos(x/10)+(x/50).^2;
%   yn = y + 0.2*randn(size(y));
%   ys = wsmooth(yn,2);
%   plot(x,yn,'.-',x,ys)
%
%   % 1-D example with non uniform steps and automatic L selection
%   x = sort(rand(1,1000)*100);
%   y = cos(x/10)+(x/50).^2;
%   yn = y + 0.2*randn(size(y));
%   [ys,L] = wsmooth(yn,x);
%   plot(x,yn,'.',x,ys,'r','Linewidth',2)
%
%   % 2-D example
%   [xp,yp] = deal(0:.02:1);
%   [x,y] = meshgrid(xp,yp);
%   f = exp(x+y) + sin((x-2*y)*3);
%   fn = f + (rand(size(f))-0.5);
%   fs = wsmooth(fn,2);
%   subplot(121), surf(xp,yp,fn), zlim([0 8])
%   subplot(122), surf(xp,yp,fs), zlim([0 8])
%
%   % 2-D example with missing data and non-uniform grid
%   xp = [0 sort(rand(1,48)) 1]; yp = [0 sort(rand(1,48)) 1];
%   [x,y] = meshgrid(xp,yp);
%   f = exp(x+y) + sin((x-2*y)*3);
%   f(round(rand(1,100)*2500)) = NaN; f(20:30,20:40) = NaN;
%   fn = f + (rand(size(f))-0.5);
%   fs = wsmooth(fn,xp,yp,2);
%   subplot(121), surf(xp,yp,fn), zlim([0 8])
%   subplot(122), surf(xp,yp,fs), zlim([0 8])
%
%   % comparison with inpaint_nans when filling NaN
%   [x,y] = meshgrid(0:.02:1);
%   [znan,z0] = deal(exp(x+y));
%   znan(10:25,20:35) = NaN;
%   znan(15:45,1:5) = NaN;
%   znan(35:37,20:45) = NaN;
%   z1 = inpaint_nans(znan);
%   z2 = wsmooth(znan,0);
%   subplot(121), surf(znan), zlim([1 8])
%   subplot(122), surf(z2), zlim([1 8]), title('wsmooth')
%
%   See also SMOOTH, INPAINT_NANS.
%
%   -- Damien Garcia -- 2008/03


%% Check input arguments
% ---
error(nargchk(1,4,nargin))

if ndims(z)~=2
    error('MATLAB:wsmooth:WrongZArray',...
        'First input must be a vector or a matrix.')
end
siz0 = size(z);

if nargin==1 % wsmooth(z)
    if isvector(z) % 1D
        z = z(:)';
        nx = length(z);
        x = 1:nx; y = 1;
    else % 2D
        [ny,nx] = size(z);
        x = 1:nx; y = 1:ny;
    end
elseif nargin==2
    if isscalar(x) && ~isvector(z) % wsmooth(z,L), 2D
        L = x;
        [ny,nx] = size(z);
        x = 1:nx; y = 1:ny;
    elseif isscalar(x) && isvector(z) % wsmooth(z,L), 1D
        L = x;
        z = z(:)';
        nx = length(z);
        x = 1:nx; y = 1;
    elseif isvector(x) && isvector(z) % wsmooth(z,x), 1D
        z = z(:)';
        y = 1;
    else
        error('MATLAB:wsmooth:WrongInputs',...
            'Wrong input arguments: type ''help wsmooth''.') 
    end
elseif nargin==3
    if isvector(z) && isvector(x) && isscalar(y) % wsmooth(z,x,L), 1D
        L = y;
        z = z(:)';
        y = 1;
    elseif ~isvector(z) && isvector(x) && isvector(y) % wsmooth(z,x,y), 2D
        % Nothing to declare for the moment
    else
        error('MATLAB:wsmooth:WrongInputs',...
            'Wrong input arguments: type ''help wsmooth''.') 
    end
else % wsmooth(z,x,y,L), 2D
    if isvector(z) || ~isvector(x) || ~isvector(y) || ~isscalar(L)
        error('MATLAB:wsmooth:WrongInputs',...
            'Wrong input arguments: type ''help wsmooth''.')   
    end
end

% Check size matching and monotonic properties
[ny,nx] = size(z);
if ~isequal(length(x),nx) || ~isequal(length(y),ny)
    error('MATLAB:wsmooth:XYLengthMismatch',...
        'Number of elements in X and/or Y is inadequate.')
elseif any(sign(diff(x))~=1) && any(sign(diff(x))~=-1)
    error('MATLAB:wsmooth:NonMonotonicX',...
        'X must be strictly monotonic.')
elseif ~isequal(y,1) && any(sign(diff(y))~=1) && any(sign(diff(y))~=-1)
    error('MATLAB:wsmooth:NonMonotonicY',...
        'Y must be strictly monotonic.')
end

% Nothing to do if dims(z)<3
if (isequal(y,1) && length(z)<3) || (~isequal(y,1) && any([nx ny]<3))
    zs = reshape(z,siz0);
    return
end

x = x(:);
y = y(:);
class0 = class(z);
z = double(z);

% Warning messages may appear in some atypical cases  
warn0 = warning('query','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:nearlySingularMatrix')

%% Create the difference matrix
%--
% 1D
if isequal(y,1)
    % Create the X-difference matrix
    % FDM: 2nd derivative, X-central difference, 2nd order
    h0 = x(2:end-1)-x(1:end-2);
    h1 = x(3:end)-x(2:end-1);    
    D = spdiags([2./h0./(h0+h1) -2./h0./h1 2./h1./(h0+h1)],...
        [0 1 2],nx-2,nx)*mean(h0.^2);
    D = D'*D;

%--    
% 2D
else
    % -----
    % The Whittaker's method is adapted to 2D. We solve:
    %   (I + lambda(Dx' Dx + Dy' Dy))zs = z,
    % where Dx and Dy are the difference matrices related to x- and y-
    % direction, respectively.
    % -----
    
    % Create the X-Difference matrix
    % FDM: 2nd derivative, X-central difference, 2nd order
    n = (nx-2)*ny;
    [Jd,Dx] = deal(zeros(3,n));

    zI = true(ny,nx);
    zI(:,[1 nx]) = false;
    [I,J] = find(zI);

    Jd(1,:) = I+(J-2)*ny;
    Jd(2,:) = I+(J-1)*ny;
    Jd(3,:) = I+J*ny;
        
    h0 = x(J)-x(J-1);
    h1 = x(J+1)-x(J);
    Dx(1,:) = 2./h0./(h0+h1);
    Dx(2,:) = -2./h0./h1;
    Dx(3,:) = 2./h1./(h0+h1);
    
    Id = repmat(1:n,[3 1]);
    D = sparse(Id,Jd,Dx,n,nx*ny)*mean(h0.^2);
    clear Dx
    D = D'*D;

    % Create the Y-Difference matrix
    % FDM: 2nd derivative, Y-central difference, 2nd order
    n = nx*(ny-2);
    [Jd,Dy] = deal(zeros(3,n));

    zI = true(ny,nx);
    zI([1 ny],:) = false;
    [I,J] = find(zI);

    Jd(1,:) = I-1+(J-1)*ny;
    Jd(2,:) = I+(J-1)*ny;
    Jd(3,:) = I+1+(J-1)*ny;

    h0 = y(I)-y(I-1);
    h1 = y(I+1)-y(I);
    Dy(1,:) = 2./h0./(h0+h1);
    Dy(2,:) = -2./h0./h1;
    Dy(3,:) = 2./h1./(h0+h1);
    
    Id = repmat(1:n,[3 1]);
    Dy = sparse(Id,Jd,Dy,n,nx*ny)*mean(h0.^2);
    clear I* J* h*
    D = D + Dy'*Dy;
    clear Dy
end

%% Weight function (0 for NaN/Inf data, 1 otherwise)
% --
% NaN or Inf data are considered missing data

I = isfinite(z);
z(~I) = 0;
W = spdiags(double(I(:)),0,nx*ny,nx*ny);

%% Solve the linear system
if exist('L','var')
    % L parameter is given
    % --------------------
    zs = (W + 10^L*D)\z(:);
else
    % Automatic determination of L parameter
    % --------------------------------------
    % The generalized cross-validation (GCV) method is used. We determine L
    % that minimizes the GCV score. If zs = H*z (with z containing N
    % elements) then the GCV score is written as:
    %        GCVs = 1/N*sum((zs-z).^2) / (1 - Tr(H)/N)^2
    % Tr(H) (trace of H) is required. Here, H = (I + L*D)^(-1). I+L*D is
    % hermitian positive semi-definite, I+L*D can be thus decomposed as
    % I+L*D = U*S*U' with U orthogonal and, S diagonal and containing the
    % eigenvalues. The eigenvalue decomposition gives: H = U*S^(-1)*U^(-1).
    % Finally: Tr(H) = Tr(inv(S)) = sum(diag(1./S)) = sum(1./eig(I+L*D)),
    % with eig(I+L*D) = 1+L*eig(D)
    
    % Note: the weight function is not taken into account in order to
    % simplify the algorithm and make it faster. This does not change the
    % results significantly.
    eigD = eig(D); % this may take a couple of seconds if D is large
    
    L = fminsearch(@(L) gcv(L),1,...
        optimset('TolX',0.1,'Display','off','Maxiter',50));
    L = round(L*10)/10;
    zs = (W + 10^L*D)\z(:);  
end

zs = reshape(zs,siz0);
zs = cast(zs,class0);
warning(warn0)

    function score = gcv(L)
    % GCV score as a function of L
    L = L(1);
    trH = sum(1./(1+10^L*eigD));
    zs = (W + 10^L*D)\z(:);
    score = sum((zs(:)-z(:)).^2)/nx/ny/(1-trH/nx/ny)^2;
    end

end

