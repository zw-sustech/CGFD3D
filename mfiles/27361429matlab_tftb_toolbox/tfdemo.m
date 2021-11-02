%TFDEMO	Demonstrate some of Time-Frequency Toolbox's capabilities.
%	TFDEMO, by itself, presents a menu of demos.

%	O.Lemoine - May 1996.
%	Copyright (c) CNRS.

clg; clear; clc;
TFTBcontinue=1;

while TFTBcontinue==1,
 choice = menu('Time-Frequency Toolbox Demonstrations',...
    'Introduction',...	
    'Non stationary signals', ...
    'Linear time-frequency representations', ...
    'Cohen''s class time-frequency distributions', ...
    'Affine class time-frequency distributions', ...
    'Reassigned time-frequency distributions', ...
    'Post-processing ',...
    'Close');
 
 if choice==1,
  tfdemo1;
 elseif choice==2,
  tfdemo2;
 elseif choice==3,
  tfdemo3;
 elseif choice==4,
  tfdemo4;
 elseif choice==5,
  tfdemo5;
 elseif choice==6,
  tfdemo6;
 elseif choice==7,
  tfdemo7;
 elseif choice==8,
  TFTBcontinue=0; 
 end
end
