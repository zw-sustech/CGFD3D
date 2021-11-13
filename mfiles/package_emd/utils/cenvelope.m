%CENVELOPE  computes envelope curves for bivariate EMD
%
% [env, mean] = CENVELOPE(t,x,ndirs,interp) computes envelope curves for bivariate EMD [1] 
%
% inputs : - x: analyzed signal
%          - t (optional): sampling times, default 1:length(x)
%          - ndirs: number of directions used to compute the mean (default: 4)
%              rem: the actual number of directions according to the paper is 2*ndirs
%          - interp (optional): interpolation scheme: 'linear', 'cubic' or 'spline' (default)
%
% outputs : - env: each stands for an envelope curve sustaining the tube envelope of
%                  the complex signal
%           - mean = mean of the envelope curves (corresponding to the first algorithm in the paper)
%
% use: env = cenvelope(x);
%      env = cenvelope(x,ndirs);
%      env = cenvelope(t,x,ndirs);
%      env = cenvelope(x,8,'linear');
%      env = cenvelope(t,x,8,'cubic');
%
%
% [1] G. Rilling, P. Flandrin, P. Gonçalves and J. M. Lilly.,
% "Bivariate Empirical Mode Decomposition",
% Signal Processing Letters (submitted)
%
% See also
%  cemd_disp (visualization)
%  emd (EMD and bivariate EMD)
%  cemdc, cemdc_fix, cemdc2, cemdc2_fix (fast implementations of bivariate EMD)
%
% G. Rilling, last modification 3.2007
% gabriel.rilling@ens-lyon.fr
function [env,envmoy] = cenvelope(t,x,Nphases,INTERP)

NBSYM = 2;
DEF_INTERP = 'spline';

if nargin < 2
    x = t;
    t = 1:length(x);
end

if nargin >= 2
    if isscalar(x)
        if nargin == 3
          INTERP = Nphases;
        end
        Nphases = x;
        x = t;
        t = 1:length(x);
    end
end

if ~exist('INTERP','var')
    INTERP = DEF_INTERP;
end
if ~exist('Nphases','var')
    Nphases = 4;
end

for k = 1:Nphases
    phi = (k-1)*pi/Nphases;
    y = real(exp(-i*phi)*x);
    [im,iM] = extr(y);
    [tmin,tmax,zmin,zmax] = boundary_conditions_emd(im,iM,t,y,x,NBSYM);
    envmin(k,:) = interp1(tmin,zmin,t,INTERP);
    envmax(k,:) = interp1(tmax,zmax,t,INTERP);
end

env=[envmin;envmax];
    
if nargout == 2
  envmoy = mean(envmax + envmin,1)/2;
end

end
