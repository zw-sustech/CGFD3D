function S=stf_rickerderiv(t,fc,t0)
%--------------------------------------------
%-- Discription:  To generate ricker wavelet
%-- Input:        t: time
%                 fc: dominant frequency
%                 t0: initial time
%-- Output:       S
%--------------------------------------------
u=(t-t0)*2.0*pi*fc;
S=u.*(-3*pi*fc+1/2*pi*fc*u.^2).*exp(-u.^2.0/4.0);
end
