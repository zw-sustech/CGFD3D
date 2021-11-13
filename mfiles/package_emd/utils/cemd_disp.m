%CEMD_DISP  displays complex envelope curves and the corresponding tube envelope
%
%
%CEMD_DISP(T,X,ENV,MODE)
%
% inputs:
%       - T: time instants
%       - X: analyzed signal (complex)
%       - ENV: matrix returned by cenvelope.m Each line is an envelope curve
%       - MODE: 'render' -> 3D rendering of the tube enclosing the signal
%               'wire'   -> wireframe display (much faster if there is no hardware acceleration)
%
% use: CEMD_DISP(X,ENV)
%      CEMD_DISP(T,X,ENV)
%      CEMD_DISP(X,ENV,MODE)
%      CEMD_DISP(T,X,ENV,MODE)
%
% rem: the program uses zbuffer rendering by default. You may switch to openGL
% rendering by commenting/uncommenting a line at the beginning of the function.
% openGL rendering is faster but zbuffer is generally nicer.
%
% See also
%  cenvelope
%
% G. Rilling, last modification: 3.2007
% gabriel.rilling@ens-lyon.fr

function varargout=cemd_disp(t,x,env,mode)

renderer = 'zbuffer'; % nicer
% renderer = 'openGL'; % faster

DEF_mode = 'render';

if nargin == 2
    env = x;
    x = t;
    t = 1:length(x);
end

if nargin == 3
    if ischar(env)
        mode = env;
        env = x;
        x = t;
        t = 1:length(x);
    end
end

if ~exist('mode','var')
    mode = DEF_mode;
end

if strcmpi(mode,'render')
    target_faces_number = 5000;
else
    target_faces_number = 1000;
end

col = [.5,.5,.5];
figure('name','complex envelope')
plot3c(t,x)
hold on
grid

[m,n] = size(env);

step = max(round(numel(env)/target_faces_number),1);
env = [env;env(1:3,:)];
inds = 1:step:n;
if inds(end)<n
    inds = [inds,n];
end
PYZ = env(:,inds);
PX = repmat(t(inds),m+3,1);

h = surf(PX,real(PYZ),imag(PYZ),col);
VN = get(h,'VertexNormals');
VN = VN(2:end-1,:,:);
PYZ = PYZ(2:end-1,:);
PX = PX(2:end-1,:);
delete(h);
h = surf(PX,real(PYZ),imag(PYZ),col,'VertexNormals',VN);
hidden off

switch lower(mode(1:min(4,end)))
    case 'rend'
        set(h,'EdgeColor','none')
        set(h,'FaceColor',col)
        for k = 1:size(env(1:end-3,:),1)
            plot3c(t,env(k,:),'k')
        end
    case 'wire'
        set(h,'EdgeColor',col)
        set(h,'FaceColor','none')
    otherwise
        warning('cemd_disp:input',['unknown option ',mode])
end
set(h,'FaceLighting','phong')
set(gcf,'renderer',renderer)

l(1) = light('Position',[min(PX(:)),min(real(PYZ(:))),min(imag(PYZ(:)))]);
l(2) = light('Position',[max(PX(:)),max(real(PYZ(:))),max(imag(PYZ(:)))]);

if nargout > 0
    varargout = {h,l};
else
    varargout = {};
end
end

function val = mod2(n,p)
val = mod(n,p);
if val == 0
    val = p;
end
end

