function S=stf_rickerderiv(t,fc,t0)

f0=sqrt(pi)/2.0;
u=(t-t0)*2.0*pi*fc;
S=u.*(1.5-u.^2/4.0).*exp(-u.^2.0/4.0)*f0*pi*fc;

