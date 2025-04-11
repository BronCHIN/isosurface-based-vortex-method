clc;
clearvars;
close all;

load('F:\work in NKRI\20170501-�������\z20171120-isosurface��ȡ��������\matlab0712\DNS_parameters.mat');
%'alpha','beta','Nx','Ny','Nz','utau','nu'; 2/3  2/3   %NX=2048 (2/3*8pi), NY=385 (1), NZ=1536 (2/3*3pi)
DetX=2/3*8*pi/Nx*utau/nu;
DetZ=2/3*3*pi/Nz*utau/nu;

load('F:\work in NKRI\20170501-�������\z20171120-isosurface��ȡ��������\matlab0712\2D_3Dlmdci_mean_max_rms.mat');
%'ind','y_plus','mean_2D_3D_lmdci','max_2D_3D_lmdci','rms_2D_3D_lmdci';
%'\lambda_{ci YZ}^{rms}','\lambda_{ci XZ}^{rms}','\lambda_{ci%XY}^{rms}','\lambda_{ci 3D}^{rms}');

%% ######### input parameters��Ҫ�޸ĵĲ���
outputnum=1000; %fprint ���ʹ��
numpolygon=60000; %Ԥ�������θ���



yic=61;%indny(kyic); %�����λ�ã����ĸ�yƽ����ȡ����
tempy=[1,5:5:190];
indtky=find(tempy<=yic);
tempky=tempy(indtky(end));
filename_lmdci=['F:\work in NKRI\20170501-�������\z20171120-isosurface��ȡ��������\data_matlab\lmd_ci_3D_Nx1536_y',num2str(tempky),'-',num2str(tempky+4),'_Nz2048.mat']
%'x_up','x_dn','y_up','y_dn','z_up','z_dn','y_step','lmdci_volm');


%% load invariables
tic;
load(filename_lmdci);
fprintf('              load data is %fs\n', toc);
fprintf('              data size=%fw\n', (x_up-x_dn+1)*(y_up-y_dn+1)*(z_up-z_dn+1)/1e4);
maxvalue=max(max_2D_3D_lmdci(y_dn:y_up,4));

for threshold=0.05:0.05:0.5 %isosurface ������ֵ
    filesavedir=['F:\work in NKRI\20170501-�������\z20171120-isosurface��ȡ��������\results\polygon_yic',num2str(yic),'_threshold',num2str(threshold),'.mat'];
    polylinez=cell(numpolygon,1);
    polylinex=cell(numpolygon,1);
    %% ����isosurface��
    tic;
    
    [Z,X,Y] = meshgrid(z_dn:z_up,x_dn:x_up,y_dn:y_step:y_up);
    [faces,vertices]=subourtine_MarchingCubes(Z,X,Y,lmdci_volm,maxvalue*threshold); %�����������
    fprintf('              isosurface calculation is %fs\n', toc);
    clear X Y Z lmdci_volm;
    
    %% Define surfaces and plot them ��������λ��
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
    %-----------���б߽��ֱ�߶Σ�
    S=Surf12;
    z1=[S.vertices(S.edges(:,1),1), S.vertices(S.edges(:,2),1)]'*DetZ; %z�����
    x1=[S.vertices(S.edges(:,1),2), S.vertices(S.edges(:,2),2)]'*DetX; %x�����
    %  figure;plot(z1,x1,'r-'); axis equal;
    %-------------------------------------------ͳ���ж��ٸ���ͬ��Ͳ�ͬ�߶�
    dot1=cat(1,z1(1,:),x1(1,:));
    dot2=cat(1,z1(2,:),x1(2,:));
    linea=cat(1,dot1,dot2)'; %�߶�����
    temp1=unique(linea,'rows');
    linearray=temp1'; %���ظ����߶�
    
    dota=cat(2,dot1,dot2)'; %������
    temp2=unique(dota,'rows');
    dotarray=temp2'; %���ظ��ĵ�
    lendot=size(dotarray,2);
    lenline=size(linearray,2);
    lendotperline=round(lendot/z_up);
    dotEdgeNo=ones(lendot,1)*100; %ͳ��ÿ�����ж��������Ǻ��Լ���ɵı�
    global POLYFLAG;%��������Ƿ��Ѿ�Ϊĳ������ڵ�
    POLYFLAG=zeros(lendot,1);
    %----------------------------------ͳ��ÿ�����ж������ߣ�ɾ������С��2���ʹ���2�ģ������ý���2���ߵĶ��㹹������Σ�
    %ע�⣡�����������������������ཻ�ģ����߲����ཻ�ģ���ɾ���ˣ�û��ͳ�ƽ�ȥ
    clear z1 x1 dot1 dot2 linea temp1 dota temp2;
    tic;
    outlen=round(lendot/outputnum); %%fprintf
    numpoly=0; %����θ���
    for k=1:lendot
        if mod(k,outlen)==0
            fprintf('             ���ݵ�%f%%�Ѵ���  ��ʱ%fs\n', round(k/lendot*10000)/100,toc);
        end
        tmz=dotarray(1,k);
        tmx=dotarray(2,k);
        
        if POLYFLAG(k)==0  %��û��������ĳ�����
            ind1=find((linearray(1,:)==tmz & linearray(2,:)==tmx)==1); %ͳ��dot1�������ж��ٸ����tmz tmxһ��
            ind2=find((linearray(3,:)==tmz & linearray(4,:)==tmx)==1); %ͳ��dot2�������ж��ٸ����tmz tmxһ��
            ind3=intersect(ind1,ind2);%����
            ind4=union(ind1,ind2);%����
            ind5=setdiff(ind4,ind3); %���,��ȥ���Լ����Լ�����߶�
            
            if  length(ind5)==2  %���������ٱ�,
                %-------������ĳһ�㿪ʼ�ıպ϶���Σ����޸ĸö���������е��POLYFLAG=1
                [closez,closex,closeflag]=subourtine_polugon(k,tmz,tmx,ind1,ind5,dotarray,linearray);
                if closeflag==1
                    numpoly=numpoly+1;
                    polylinez{numpoly}=closez;
                    polylinex{numpoly}=closex;
                end
            else
                POLYFLAG(k)=2; %�����������ٱߵģ�Ҳ�ų�����
            end%if�Ƿ�������
            
        end%if�Ƿ��Ѿ������if
    end%k
    fprintf('              �����պ϶���Σ�����ʱ�� %fs\n', toc);
    
    
    %% Calculate radius circularity elongation #2
    %����Sommer ���ڵĳ�����㿴��
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

