function M=src_mech(strike,dip,rake)
dip_pi=dip/180*pi;
strike_pi=strike/180*pi;
rake_pi=rake/180*pi;

Mxx=-(sin(dip_pi)*cos(rake_pi)*sin(2.0*strike_pi)   ...
               +sin(2.0*dip_pi)*sin(rake_pi)*sin(strike_pi)^2);
Myy=sin(dip_pi)*cos(rake_pi)*sin(2.0*strike_pi)   ...
               -sin(2.0*dip_pi)*sin(rake_pi)*cos(strike_pi)^2;
Mzz=-(Mxx+Myy);
Mxy=sin(dip_pi)*cos(rake_pi)*cos(2.0*strike_pi)     ...
               +0.5*sin(2.0*dip_pi)*sin(rake_pi)*sin(2.0*strike_pi);
Mxz=-(cos(dip_pi)*cos(rake_pi)*cos(strike_pi)       ...
               +cos(2.0*dip_pi)*sin(rake_pi)*sin(strike_pi));
Myz=-(cos(dip_pi)*cos(rake_pi)*sin(strike_pi)       ...
               -cos(2.0*dip_pi)*sin(rake_pi)*cos(strike_pi));
%Mxz=-Mxz;Mxy=-Mxy !for upward positive z axis
% convert to x:east y:norht z:upward system
%Mtm=Mxx; Mxx=Myy; Myy=Mtm
%Mtm=Mxz; Mxz=-Myz; Myz=-Mtm
% convert to theta:south phi:east r:upward system
%M=[ Myy  -Mxy  -Myz;
%   -Mxy   Mxx   Mxz;
%   -Myz   Mxz   Mzz];

% convert to theta:south phi:east r:upward system
M=[ Mxx  -Mxy   Mxz;
   -Mxy   Myy  -Myz;
    Mxz  -Myz   Mzz];

