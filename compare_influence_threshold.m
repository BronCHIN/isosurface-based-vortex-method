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
for threshold=0.05:0.05:0.3 %isosurface 定义阈值
    filedir=['F:\work in NKRI\20170501-青年基金\z20171120-isosurface提取涡旋特征\matlab0713_threshold\results\polygon_yic',num2str(yic),'_threshold',num2str(threshold),'.mat'];
    load(filedir);%'numpoly','polylinez','polylinex','areapys','radiuspys','circularity','elongation');
    num=num+1;
    Novort(num)=numpoly;
    Noarea(num,1)=mean(areapys);
    Noarea(num,2)=std(areapys,1);
    Noradis(num,1)=mean(radiuspys);
    Noradis(num,2)=std(radiuspys,1);
    Nocir(num,1)=mean(circularity);
    Nocir(num,2)=std(circularity,1);
    Noenlg(num,1)=mean(elongation);
    Noenlg(num,2)=std(elongation,1);
    
    
    figure;
    for k=1:numpoly
        patch(polylinez{k},polylinex{k},'red');            %just for displaying purposes
    end
    axis equal
    
end
figure(1);
semilogy(Nothrd,Novort,'k.-')  ;
xlabel('阈值');ylabel('提取涡旋个数');hold on;

figure(2);
%errorbar(Nothrd,Noarea(:,1),Noarea(:,2),'k-o')
plot(Nothrd,Noarea(:,1),'k.-')  ;%
xlabel('阈值');ylabel('面积');hold on;

figure(3);
plot(Nothrd,Noradis(:,1),'k.-')  ;%Noradis(:,2),'k-o')
xlabel('阈值');ylabel('半径');hold on;

figure(4);
plot(Nothrd,Nocir(:,1),'k.-')  ;%Nocir(:,2),'k-o')
xlabel('阈值');ylabel('圆度');hold on;

figure(5);
plot(Nothrd,Noenlg(:,1),'k.-')  ;%Noenlg(:,2),'k-o')
xlabel('阈值');ylabel('狭长度');hold on;
% step=0.01;
% nhist=0:step:1;
% [nt1,xout1] =hist(circularity,nhist); nt1=nt1/numpoly;
% disp(['均值圆度=',num2str(mean(circularity))]);
%
% step=0.25;
% nhist=0:step:30;
% [nt2,xout2] =hist(elongation,nhist); nt2=nt2/numpoly;
% disp(['均值狭长度=',num2str(mean(elongation))]);
%
% step=50;
% nhist=0:step:2000;
% [nt3,xout3] =hist(areapys,nhist); nt3=nt3/numpoly;
% disp(['均值面积=',num2str(mean(areapys))]);
%
% step=0.1;
% nhist=0:step:40;
% [nt4,xout4] =hist(radiuspys,nhist); nt4=nt4/numpoly;
% disp(['均值半径=',num2str(mean(radiuspys))]);
%
% figure(1);plot(xout1,nt1,'k.-'); xlabel('circularity');ylabel('PDF');hold on;
% figure(2);plot(xout2,nt2,'k.-'); xlabel('elongation');ylabel('PDF');hold on;
% figure(3);plot(xout3,nt3,'k.-'); xlabel('area');ylabel('PDF');hold on;
% figure(4);plot(xout4,nt4,'k.-'); xlabel('radius');ylabel('PDF');hold on;

