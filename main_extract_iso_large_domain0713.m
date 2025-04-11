clc;
clearvars;
close all;

load('F:\work in NKRI\20170501-青年基金\z20171120-isosurface提取涡旋特征\matlab0712\DNS_parameters.mat');
%'alpha','beta','Nx','Ny','Nz','utau','nu'; 2/3  2/3   %NX=2048 (2/3*8pi), NY=385 (1), NZ=1536 (2/3*3pi)
DetX=2/3*8*pi/Nx*utau/nu;
DetZ=2/3*3*pi/Nz*utau/nu;

load('F:\work in NKRI\20170501-青年基金\z20171120-isosurface提取涡旋特征\matlab0712\2D_3Dlmdci_mean_max_rms.mat');
%'ind','y_plus','mean_2D_3D_lmdci','max_2D_3D_lmdci','rms_2D_3D_lmdci';
%'\lambda_{ci YZ}^{rms}','\lambda_{ci XZ}^{rms}','\lambda_{ci%XY}^{rms}','\lambda_{ci 3D}^{rms}');

%% ######### input parameters需要修改的参数
outputnum=1000; %fprint 语句使用
numpolygon=60000; %预定义多边形个数



yic=61;%indny(kyic); %切面的位置，从哪个y平面提取涡旋
tempy=[1,5:5:190];
indtky=find(tempy<=yic);
tempky=tempy(indtky(end));
filename_lmdci=['F:\work in NKRI\20170501-青年基金\z20171120-isosurface提取涡旋特征\data_matlab\lmd_ci_3D_Nx1536_y',num2str(tempky),'-',num2str(tempky+4),'_Nz2048.mat']
%'x_up','x_dn','y_up','y_dn','z_up','z_dn','y_step','lmdci_volm');


%% load invariables
tic;
load(filename_lmdci);
fprintf('              load data is %fs\n', toc);
fprintf('              data size=%fw\n', (x_up-x_dn+1)*(y_up-y_dn+1)*(z_up-z_dn+1)/1e4);
maxvalue=max(max_2D_3D_lmdci(y_dn:y_up,4));

for threshold=0.05:0.05:0.5 %isosurface 定义阈值
    filesavedir=['F:\work in NKRI\20170501-青年基金\z20171120-isosurface提取涡旋特征\results\polygon_yic',num2str(yic),'_threshold',num2str(threshold),'.mat'];
    polylinez=cell(numpolygon,1);
    polylinex=cell(numpolygon,1);
    %% 计算isosurface面
    tic;
    
    [Z,X,Y] = meshgrid(z_dn:z_up,x_dn:x_up,y_dn:y_step:y_up);
    [faces,vertices]=subourtine_MarchingCubes(Z,X,Y,lmdci_volm,maxvalue*threshold); %算起来超快的
    fprintf('              isosurface calculation is %fs\n', toc);
    clear X Y Z lmdci_volm;
    
    %% Define surfaces and plot them 定义切面位置
    % Create Surface #1
    Surface1.vertices = [z_dn,x_dn-round((x_up+x_dn)/2),yic; z_dn,x_up+round((x_up+x_dn)/2),yic; z_up+(z_up-z_dn),round((x_up+x_dn)/2),yic];
    Surface1.faces    = [1 2 3];
    % Create Surface #2
    Surface2.vertices=vertices;
    Surface2.faces=faces;
    
    %% Run SurfaceIntersection
    tic;
    [~, Surf12] = subourtine_SurfaceIntersection(Surface1, Surface2);
    fprintf('              intersection calculation is %fs\n', toc);
    
    %% Extract vortex
    %-----------所有边界的直线段，
    S=Surf12;
    z1=[S.vertices(S.edges(:,1),1), S.vertices(S.edges(:,2),1)]'*DetZ; %z方向的
    x1=[S.vertices(S.edges(:,1),2), S.vertices(S.edges(:,2),2)]'*DetX; %x方向的
    %  figure;plot(z1,x1,'r-'); axis equal;
    %-------------------------------------------统计有多少个不同点和不同线段
    dot1=cat(1,z1(1,:),x1(1,:));
    dot2=cat(1,z1(2,:),x1(2,:));
    linea=cat(1,dot1,dot2)'; %线段序列
    temp1=unique(linea,'rows');
    linearray=temp1'; %不重复的线段
    
    dota=cat(2,dot1,dot2)'; %点序列
    temp2=unique(dota,'rows');
    dotarray=temp2'; %不重复的点
    lendot=size(dotarray,2);
    lenline=size(linearray,2);
    lendotperline=round(lendot/z_up);
    dotEdgeNo=ones(lendot,1)*100; %统计每个点有多少条不是和自己组成的边
    global POLYFLAG;%用来标记是否已经为某多边形内点
    POLYFLAG=zeros(lendot,1);
    %----------------------------------统计每个点有多少条边（删除边数小于2，和大于2的），并用仅有2条边的顶点构建多边形，
    %注意！！！！！！！！这里面有相交的，或者部分相交的，就删掉了，没有统计进去
    clear z1 x1 dot1 dot2 linea temp1 dota temp2;
    tic;
    outlen=round(lendot/outputnum); %%fprintf
    numpoly=0; %多边形个数
    for k=1:lendot
        if mod(k,outlen)==0
            fprintf('             数据第%f%%已处理  耗时%fs\n', round(k/lendot*10000)/100,toc);
        end
        tmz=dotarray(1,k);
        tmx=dotarray(2,k);
        
        if POLYFLAG(k)==0  %且没被归类于某多边形
            ind1=find((linearray(1,:)==tmz & linearray(2,:)==tmx)==1); %统计dot1矩阵里有多少个点和tmz tmx一致
            ind2=find((linearray(3,:)==tmz & linearray(4,:)==tmx)==1); %统计dot2矩阵里有多少个点和tmz tmx一致
            ind3=intersect(ind1,ind2);%并集
            ind4=union(ind1,ind2);%交集
            ind5=setdiff(ind4,ind3); %差集合,除去点自己和自己组成线段
            
            if  length(ind5)==2  %仅有两条临边,
                %-------计算以某一点开始的闭合多边形，并修改该多边形上所有点的POLYFLAG=1
                [closez,closex,closeflag]=subourtine_polugon(k,tmz,tmx,ind1,ind5,dotarray,linearray);
                if closeflag==1
                    numpoly=numpoly+1;
                    polylinez{numpoly}=closez;
                    polylinex{numpoly}=closex;
                end
            else
                POLYFLAG(k)=2; %不符合两条临边的，也排除在外
            end%if是否两条边
            
        end%if是否已经多边形if
    end%k
    fprintf('              构建闭合多边形，花费时间 %fs\n', toc);
    
    
    %% Calculate radius circularity elongation #2
    %利用Sommer 教授的程序计算看看
    areapys=zeros(numpoly,1);
    perimeter=zeros(numpoly,1);
    circularity=zeros(numpoly,1);
    elongation=zeros(numpoly,1);
    
    for k=1:numpoly
        [ geom, iner, cpmo ] = subourtine_polygeom(polylinez{k},polylinex{k}) ;
        % GEOM = [ area   X_cen  Y_cen  perimeter ]
        %INER = [ Ixx    Iyy    Ixy    Iuu    Ivv    Iuv ]
        %CPMO = [ I1     ang1   I2     ang2   J ]
        areapys(k) = geom(1);
        circularity(k)=4*pi*geom(1)./geom(4)/geom(4);
        elongation(k)=sqrt(cpmo(3)/cpmo(1));
    end
    radiuspys=sqrt(areapys/pi);
    
    %% output
    save(filesavedir,'numpoly','polylinez','polylinex','areapys','radiuspys','circularity','elongation');
    clear faces vertices polylinez  polylinex  linearray dotarray areapys radiuspys circularity elongation;
end

