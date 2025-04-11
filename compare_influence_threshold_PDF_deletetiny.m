%对比计算区域对结果的影响
clc;
clear all;
close all;
load('F:\work in NKRI\20170501-青年基金\z20171120-isosurface提取涡旋特征\matlab0712\DNS_parameters.mat');
%'alpha','beta','Nx','Ny','Nz','utau','nu'; 2/3  2/3   %NX=2048 (2/3*8pi), NY=385 (1), NZ=1536 (2/3*3pi)
%三维空间坐标x,y,z
yp=193;
x = (0:(Nx-1))/(Nx)*(8.0*pi*alpha);
y = 1-cos(pi*(0:(yp-1))/(Ny-1));
z = (0:(Nz-1))/(Nz)*(3.0*pi*beta);
y_plus=y*utau/nu;
x_plus=x*utau/nu;

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

num=0;
nummark=1;
clor='kbgrcym';
mark='.+xsdo^';
for threshold=0.05:0.05:0.5 %isosurface 定义阈值
    filedir=['F:\work in NKRI\20170501-青年基金\z20171120-isosurface提取涡旋特征\matlab0713_threshold\results\polygon_yic',num2str(yic),'_threshold',num2str(threshold),'_delete_tiny.mat'];
    load(filedir);%'numpoly2','polylinez','polylinex','areapys','radiuspys','circularity','elongation');
    num=num+1;
    if num>7
        nummark=nummark+1;
    end
    Novort(num)=numpoly2;
    Noarea(num,1)=mean(areapys2);
    Noarea(num,2)=std(areapys2,1);
    Noradis(num,1)=mean(radiuspys2);
    Noradis(num,2)=std(radiuspys2,1);
    Nocir(num,1)=mean(circularity2);
    Nocir(num,2)=std(circularity2,1);
    Noenlg(num,1)=mean(elongation2);
    Noenlg(num,2)=std(elongation2,1);
    
    step=0.05;
    nhist=0:step:1;
    [nt1,xout1] =hist(circularity2,nhist); nt1=nt1/numpoly2;
    figure(1);plot(xout1,nt1,[clor(min(num,7)),mark(nummark),'-']); xlabel('circularity');ylabel('PDF');hold on;
    
    step=0.25;
    nhist=0:step:30;
    [nt2,xout2] =hist(elongation2,nhist); nt2=nt2/numpoly2;
     figure(2);plot(xout2,nt2,[clor(min(num,7)),mark(nummark),'-']); xlabel('elongation');ylabel('PDF');hold on;
     
    step=50;
    nhist=0:step:2000;
    [nt3,xout3] =hist(areapys2,nhist); nt3=nt3/numpoly2;
    figure(3);plot(xout3,nt3,[clor(min(num,7)),mark(nummark),'-']); xlabel('area');ylabel('PDF');hold on;
    
    step=0.5;
    nhist=0:step:40;
    [nt4,xout4] =hist(radiuspys2,nhist); nt4=nt4/numpoly2;
    figure(4);plot(xout4,nt4,[clor(min(num,7)),mark(nummark),'-']); xlabel('radius');ylabel('PDF');hold on;
end
figure(1);legend('0.05','0.1','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50');
figure(2);legend('0.05','0.1','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50');
figure(3);legend('0.05','0.1','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50');
figure(4);legend('0.05','0.1','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50');

%
% figure(1);plot(xout1,nt1,'k.-'); xlabel('circularity');ylabel('PDF');hold on;
% figure(2);plot(xout2,nt2,'k.-'); xlabel('elongation');ylabel('PDF');hold on;
% figure(3);plot(xout3,nt3,'k.-'); xlabel('area');ylabel('PDF');hold on;
% figure(4);plot(xout4,nt4,'k.-'); xlabel('radius');ylabel('PDF');hold on;

