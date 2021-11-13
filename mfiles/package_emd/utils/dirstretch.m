%DIRSTRETCH  variable directional stretching of a complex vector
%
% DIRSTRETCH(X,D) where X is a complex vector and D a complex scalar stretches X by abs(D) along the direction arg(D)
%   if D is a vector the same size as X, then X(i) is stretched by abs(D(i)) along the direction arf(D(i))
% DIRSTRETCH(X,DIR,COEF) where X is a complex vector, DIR is either a unit modulus complex number or a real number and COEF a real number 
%   stretches X  by COEF along the direction DIR (or EXP(I*DIR) if DIR is real)
%   if DIR and/or COEF is a vector the same size as X, then X(i) is stretched by COEF(i) along the direction DIR(i) (or EXP(I*DIR(i)))
%
% G. Rilling 3.2007
% gabriel.rilling@ens-lyon.fr
function y=dirstretch(x,dir,coef)
if nargin == 3 
  if ~isreal(dir) && all(abs(dir) ~= 1)
    warning('non unit modulus of direction argument ignored')
  end
  if isreal(dir) && any(abs(dir) ~= 1)
    dir = exp(i*dir);
  end
  if any(dir==0)
    error('invalid zero direction argument')
  end
  tmp = abs(dir);
  tmp(tmp==0) = 1; % to avoid NaNs
  dir = dir./tmp;
end
if nargin == 2
    coef = abs(dir);
    tmp = coef;
    tmp(coef==0) = 1; % to avoid NaNs
    dir = 1./tmp.*dir;
end
rotx = x.*conj(dir);
y = coef.*dir.*real(rotx)+i*dir.*imag(rotx);