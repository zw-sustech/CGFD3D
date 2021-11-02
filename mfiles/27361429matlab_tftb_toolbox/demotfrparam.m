% programme de test

NbPoints=600; 
NbPoints-rem(NbPoints,2); % force NbPoints to be even

bruit=randn(NbPoints,1);
signal1=filter([1 0 0],[1 -2*0.98*cos(2*pi*0.1) 1],bruit);
signal1=signal1(NbPoints/2+1:NbPoints); signal1=signal1/std(signal1);
signal2=filter([1 0 0],[1 -2*0.98*cos(2*pi*0.3) 1],bruit);
signal2=signal2(NbPoints/2+1:NbPoints); signal2=signal2/std(signal2);
signal=[signal1;signal2];

figure(1); plot(signal); drawnow;

[psd1,freqs]=parafrep(correlmx(signal( 10+(1:50)),2,'fbhermitian'),128,'ar');
[psd2,freqs]=parafrep(correlmx(signal(500+(1:50)),2,'fbhermitian'),128,'ar');

figure(2); plot(freqs,20*log10(psd1),freqs,20*log10(psd2)); drawnow;

figure(3); tfrparam(signal,10:5:600,256,2,51,1,'fbhermitian','ar',1);
figure(4); tfrparam(signal,10:5:600,256,2,51,1,'fbhermitian','capon',1);
figure(5); tfrparam(signal,10:5:600,256,2,51,1,'fbhermitian','lagunas',1);
figure(6); tfrparam(signal,10:5:600,256,30,51,1,'fbhermitian','periodogram',1);


