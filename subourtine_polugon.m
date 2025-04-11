function [closez,closex,closeflag]=subourtine_polugon(k,tmz,tmx,ind1,ind5,dotarray,linearray)
%written by Bron 20180709-10
global POLYFLAG; %用来标记是否已经为某多边形内点
POLYFLAG(k,1)=1; %两条临边
closez=[]; closez=[closez tmz];
closex=[]; closex=[closex tmx];
closeflag=0; %初始化，不闭合

tempzx=linearray(:,ind5); %计算两条临边的其他两点坐标
if ismember(ind5(1),ind1) %看ind5中第一个元素 在 ind1中
    startz=tempzx(3,1); %ind1->(3,4)
    startx=tempzx(4,1);%ind1
else  %看ind5中第一个元素 在 ind2中
    startz=tempzx(1,1);%ind2->(1,2)
    startx=tempzx(2,1);
end
if ismember(ind5(2),ind1) %看ind5中第二个元素 在 ind1中
    endz0=tempzx(3,2);%ind1->(3,4)
    endx0=tempzx(4,2);
else  %看ind5中第一个元素 在 ind2中
    endz0=tempzx(1,2);%ind2->(1,2)
    endx0=tempzx(2,2);
end


while (startz~=endz0 && startx~=endx0)
    
    ind1=find((linearray(1,:)==startz & linearray(2,:)==startx)==1); %统计dot1矩阵里有多少个点和tmz tmx一致
    ind2=find((linearray(3,:)==startz & linearray(4,:)==startx)==1); %统计dot2矩阵里有多少个点和tmz tmx一致
    ind3=intersect(ind1,ind2);%并集
    ind4=union(ind1,ind2);%交集
    ind5=setdiff(ind4,ind3); %差集合,除去点自己和自己组成线段
    
    if  length(ind5)==2  %仅有两条临边,
        %--封闭多边形沿线坐标
        closez=[closez startz];
        closex=[closex startx];
        %----修改点属性多边形范围
        ind6=(dotarray(1,:)==startz & dotarray(2,:)==startx);
        POLYFLAG(ind6)=1;
        
        tempzx=linearray(:,ind5); %计算两条临边的其他两点坐标
        %-----------计算点的另外相邻两点
        if ismember(ind5(1),ind1) %看ind5中第一个元素 在 ind1中
            tz1=tempzx(3,1); %ind1->(3,4)
            tx1=tempzx(4,1);%ind1
        else  %看ind5中第一个元素 在 ind2中
            tz1=tempzx(1,1);%ind2->(1,2)
            tx1=tempzx(2,1);
        end
        if ismember(ind5(2),ind1) %看ind5中第二个元素 在 ind1中
            tz2=tempzx(3,2);%ind1->(3,4)
            tx2=tempzx(4,2);
        else  %看ind5中第一个元素 在 ind2中
            tz2=tempzx(1,2);%ind2->(1,2)
            tx2=tempzx(2,2);
        end
        %------闭合路线方向选择
        if tz1==tmz && tx1==tmx %看哪个点是之前的点，
            nextz=tz2;
            nextx=tx2;
        elseif tz2==tmz && tx2==tmx
            nextz=tz1;
            nextx=tx1;
        end
        %----前一个点
        tmz=startz;
        tmx=startx;
        %---下一个点
        startz=nextz;
        startx=nextx;
        
    else
        ind6=(dotarray(1,:)==startz & dotarray(2,:)==startx);
        POLYFLAG(ind6)=2; %不符合两条临边的，也排除在外
        return;%直接结束程序
    end % 是否两条边
end %while

if (startz==endz0 && startx==endx0) && length(ind5)==2 %曲线封闭
    %--封闭多边形最后一点坐标
    closez=[closez endz0];
    closex=[closex endx0];
    %----修改点属性多边形范围
    ind6=(dotarray(1,:)==endz0 & dotarray(2,:)==endx0);
    POLYFLAG(ind6)=1;
    
    closeflag=1;
end

end