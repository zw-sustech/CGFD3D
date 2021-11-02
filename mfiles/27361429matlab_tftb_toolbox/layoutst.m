

issig=0;
isspec=0;
iscolorbar=0;

linlogspec=1;
sigenveloppe=0;


layout=issig+isspec*2+iscolorbar*4;
while layout~=4,
 if issig==0, 
  SignalStr= 'display signal';
 else
  SignalStr='remove signal';
 end;
 
 if isspec==0, 
  SpectrumStr= 'display spectrum';
 else
  SpectrumStr='remove spectrum';
 end;
 
 if iscolorbar==0,
  ColorbarStr='display colorbar';
 else
  ColorbarStr='remove colorbar';
 end;

 layout=menu('DISPLAY LAYOUT',...
             SignalStr,...
             SpectrumStr,...
             ColorbarStr,...
             'close');
            
 if layout==1,
  issig=~issig;
  if issig==1,
   sigenveloppe=menu('SIGNAL REPRESENTATION','signal only','signal with enveloppe')-1;
  end; 
 elseif layout==2,
  isspec=~isspec;
  if isspec==1,
   linlogspec=menu('FREQUENCY REPRESENTATION','linear scale','log scale')-1;
  end;
 elseif layout==3,
  iscolorbar=~iscolorbar;
 end;             
end;           




