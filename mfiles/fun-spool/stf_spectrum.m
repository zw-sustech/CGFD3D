function [f,F]=stf_spectrum(y,Fs);

L=length(y);

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2);

F=2*abs(Y(1:NFFT/2));
