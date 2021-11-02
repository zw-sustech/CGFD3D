% EX_ONLINE.M
%
% P. Flandrin, Mar. 13, 2003
%
% computes on-line EMD for the periodized version of the
% sum of 2 sinusoidal FM's + 1 Gaussian logon used in EMD_FMSIN.M
%
% displays reassigned spectrograms of the sum signal and of the 3 first
% modes extracted
%
% illustrates Section 3.4 in
%
% G. Rilling, P. Flandrin and P. Gonçalvès
% "On Empirical Mode Decomposition and its algorithms"
% IEEE-EURASIP Workshop on Nonlinear Signal and Image Processing
% NSIP-03, Grado (I), June 2003

% elementary signal
N = 2000;% # of data samples
T = 1:4:N;
t = 1:N;

p = N/2;% period of the 2 sinusoidal FM's

% sinusoidal FM 1
fmin1 = 1/64;% min frequency
fmax1 = 1.5*1/8;% max frequency
x1 = fmsin(N,fmin1,fmax1,p,N/2,fmax1);

% sinusoidal FM 2
fmin2 = 1/32;% min frequency
fmax2 = 1.5*1/4;% max frequency
x2 = fmsin(N,fmin2,fmax2,p,N/2,fmax2);

% logon
f0 = 1.5*1/16;% center frequency
x3 = amgauss(N,N/2,N/8).*fmconst(N,f0);

a1 = 1;
a2 = 1;
a3 = 1;

x = real(a1*x1+a2*x2+a3*x3);
x = x/max(abs(x));

% periodized signal
y = [x' x' x' x' x' x' x' x'];

[imf_el,ort_el,nbits_el] = emd_online(y,1:length(y),[0.05,0.5,20,100],4,1,1);

T = 1:32:length(y);

emd_visu(y,1:length(y),imf_el,1);

figure(1)

% time-frequency distributions
Nf = 256;% # of frequency bins
Nh = 127;% short-time window length
w = window(Nh,'Kaiser');

[s,rs] = tfrrsp(y',T,Nf,w,1);
[s,rs1] = tfrrsp(imf_el(1,:)',T,Nf,w,1);
[s,rs2] = tfrrsp(imf_el(2,:)',T,Nf,w,1);
[s,rs3] = tfrrsp(imf_el(3,:)',T,Nf,w,1);

figure(4)

subplot(221)
imagesc(flipud(rs(1:128,:)))
set(gca,'YTick',[]);set(gca,'XTick',[])
xlabel('time')
ylabel('frequency')
title('signal')

subplot(222)
imagesc(flipud(rs1(1:128,:)))
set(gca,'YTick',[]);set(gca,'XTick',[])
xlabel('time')
ylabel('frequency')
title('mode #1')

subplot(223)
imagesc(flipud(rs2(1:128,:)))
set(gca,'YTick',[]);set(gca,'XTick',[])
xlabel('time')
ylabel('frequency')
title('mode #2')

subplot(224)
imagesc(flipud(rs3(1:128,:)))
set(gca,'YTick',[]);set(gca,'XTick',[])
xlabel('time')
ylabel('frequency')
title('mode #3')

%colormap(flipud(gray))
colormap(jet)
