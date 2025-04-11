clc;
clear all;
close all;
%%
load('F:\work in NKRI\20170501-青年基金\z20171120-isosurface提取涡旋特征\matlab0712\DNS_parameters.mat');
%'alpha','beta','Nx','Ny','Nz','utau','nu'; 2/3  2/3   %NX=2048 (2/3*8pi), NY=385 (1), NZ=1536 (2/3*3pi)
%三维空间坐标x,y,z
yp=193;
x = (0:(Nx-1))/(Nx)*(8.0*pi*alpha);
y = 1-cos(pi*(0:(yp-1))/(Ny-1));
z = (0:(Nz-1))/(Nz)*(3.0*pi*beta);
y_plus=y*utau/nu;
x_plus=x*utau/nu;
DetX=2/3*8*pi/Nx*utau/nu;
DetZ=2/3*3*pi/Nz*utau/nu;

load('F:\work in NKRI\20170501-青年基金\z20171120-isosurface提取涡旋特征\matlab0712\2D_3Dlmdci_mean_max_rms.mat');
%'ind','y_plus','mean_2D_3D_lmdci','max_2D_3D_lmdci','rms_2D_3D_lmdci';
%'\lambda_{ci YZ}^{rms}','\lambda_{ci XZ}^{rms}','\lambda_{ci%XY}^{rms}','\lambda_{ci 3D}^{rms}');

yic=61;%indny(kyic); %切面的位置，从哪个y平面提取涡旋
Nothrd=0.05:0.05:0.5;
Novort=zeros(size(Nothrd));
Noarea=zeros(length(Nothrd),2);
Noradis=zeros(length(Nothrd),2);
Noenlg=zeros(length(Nothrd),2);
Nocir=zeros(length(Nothrd),2);


for threshold=0.05:0.05:0.5 %isosurface 定义阈值
    filedir=['F:\work in NKRI\20170501-青年基金\z20171120-isosurface提取涡旋特征\matlab0713_threshold\results\polygon_yic',num2str(yic),'_threshold',num2str(threshold),'.mat'];
    load(filedir);%'numpoly','polylinez','polylinex','areapys','radiuspys','circularity','elongation');
    threshold
    %% Calculate radius circularity elongation #2
    %利用Sommer 教授的程序计算看看
   
    num=0;
    for k=1:numpoly
        difz=max(polylinez{k})-min(polylinez{k});
        difx=max(polylinex{k})-min(polylinex{k});
        flag=(difz>0.25*DetZ) && (difx>0.25*DetX);
        if flag==1
            num=num+1;
            [ geom, iner, cpmo ] = subourtine_polygeom(polylinez{k},polylinex{k}) ;
            % GEOM = [ area   X_cen  Y_cen  perimeter ]
            %INER = [ Ixx    Iyy    Ixy    Iuu    Ivv    Iuv ]
            %CPMO = [ I1     ang1   I2     ang2   J ]
            areapys2(num) = geom(1); %已经无量纲化了
            circularity2(num)=4*pi*geom(1)./geom(4)/geom(4);
            elongation2(num)=sqrt(cpmo(3)/cpmo(1));
        end
    end
    radiuspys2=sqrt(areapys2/pi);
        
    %% output
    numpoly2=num;
    filesavedir=['F:\work in NKRI\20170501-青年基金\z20171120-isosurface提取涡旋特征\matlab0713_threshold\results\polygon_yic',num2str(yic),'_threshold',num2str(threshold),'_delete_tiny.mat'];
    save(filesavedir,'numpoly','numpoly2','polylinez','polylinex','areapys2','radiuspys2','circularity2','elongation2');
    clear areapys2 radiuspys2 circularity2 elongation2;
end