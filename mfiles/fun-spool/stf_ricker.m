function S=stf_ricker(t,fc,t0)

nt=length(t);

f0=sqrt(pi)/2.0;
u=(t-t0)*2.0*pi*fc;
S=(u.^2.0/4.0-0.5).*exp(-u.^2.0/4.0)*f0;

