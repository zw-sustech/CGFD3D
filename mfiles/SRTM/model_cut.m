clc
clear

XN=990;
YN=1376;
ZN=70;


y_min=21;
x_min=97.5;


x1=103;
x2=106;
y1=28;
y2=31;

s_sy=28.34;
s_sx=104.90;
ns_sx=(s_sx-x1)*110000
ns_sy=(s_sy-y1)*110000



xn=(x2-x1)*110;
yn=(y2-y1)*110;
zn=30;

nvp=zeros(xn,yn,zn);
nvs=zeros(xn,yn,zn);
np=zeros(xn,yn,zn);

vp=zeros(XN,YN,ZN);
vs=zeros(XN,YN,ZN);
dp=zeros(XN,YN,ZN);

fp=fopen('vp.bin','rb');
vps=fread(fp,inf,'float');
fclose(fp);

fp=fopen('vs.bin','rb');
vss=fread(fp,inf,'float');
fclose(fp);

fp=fopen('p.bin','rb');
ps=fread(fp,inf,'float');
fclose(fp);
 for k=1:ZN
    for j=1:YN
        for i=1:XN   
            vp(i,j,k)=vps((k-1)*YN*XN+(j-1)*XN+i);
            vs(i,j,k)=vss((k-1)*YN*XN+(j-1)*XN+i);
            dp(i,j,k)=ps((k-1)*YN*XN+(j-1)*XN+i);
        end
    end
 end

 for k=1:zn
    for j=1:yn
        for i=1:xn  
            nvp(i,j,k)=vp(i+(x1-x_min)*110,j+(y1-y_min)*110,k);
            nvs(i,j,k)=vs(i+(x1-x_min)*110,j+(y1-y_min)*110,k);
            np(i,j,k)=dp(i+(x1-x_min)*110,j+(y1-y_min)*110,k);
            
        end
    end
 end
 
 

fp=fopen('YBvp.bin','wb');
fwrite(fp,nvp,'float');
fclose(fp);

fp=fopen('YBvs.bin','wb');
fwrite(fp,nvs,'float');
fclose(fp);

fp=fopen('YBp.bin','wb');
fwrite(fp,np,'float');
fclose(fp);



[xx3,yy3,zz3] = meshgrid( 0:1:(yn-1),0:1:(xn-1)*1, 0:1:(zn-1));

figure(1)
slice(xx3,yy3,zz3,np,[ 250 350],[ 113 200],[20 50]); 
shading interp
colorbar
colormap('jet')
set(gca,'zdir','reverse');
% axis equal
grid on
box on