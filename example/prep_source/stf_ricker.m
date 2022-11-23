function S=stf_ricker(t,fc,t0)
%--------------------------------------------
%-- Discription:  To generate ricker wavelet
%-- Input:        t: time
%                 fc: dominant frequency
%                 t0: offset time
%-- Output:       S
%-------------------------------------------
u = (t-t0)*2.0*pi*fc;
S = (1-u.^2.0/2.0).*exp(-u.^2.0/4.0);
end

