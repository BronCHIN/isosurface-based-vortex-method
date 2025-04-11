function [closez,closex,closeflag]=subourtine_polugon(k,tmz,tmx,ind1,ind5,dotarray,linearray)
%written by Bron 20180709-10
global POLYFLAG; %��������Ƿ��Ѿ�Ϊĳ������ڵ�
POLYFLAG(k,1)=1; %�����ٱ�
closez=[]; closez=[closez tmz];
closex=[]; closex=[closex tmx];
closeflag=0; %��ʼ�������պ�

tempzx=linearray(:,ind5); %���������ٱߵ�������������
if ismember(ind5(1),ind1) %��ind5�е�һ��Ԫ�� �� ind1��
    startz=tempzx(3,1); %ind1->(3,4)
    startx=tempzx(4,1);%ind1
else  %��ind5�е�һ��Ԫ�� �� ind2��
    startz=tempzx(1,1);%ind2->(1,2)
    startx=tempzx(2,1);
end
if ismember(ind5(2),ind1) %��ind5�еڶ���Ԫ�� �� ind1��
    endz0=tempzx(3,2);%ind1->(3,4)
    endx0=tempzx(4,2);
else  %��ind5�е�һ��Ԫ�� �� ind2��
    endz0=tempzx(1,2);%ind2->(1,2)
    endx0=tempzx(2,2);
end


while (startz~=endz0 && startx~=endx0)
    
    ind1=find((linearray(1,:)==startz & linearray(2,:)==startx)==1); %ͳ��dot1�������ж��ٸ����tmz tmxһ��
    ind2=find((linearray(3,:)==startz & linearray(4,:)==startx)==1); %ͳ��dot2�������ж��ٸ����tmz tmxһ��
    ind3=intersect(ind1,ind2);%����
    ind4=union(ind1,ind2);%����
    ind5=setdiff(ind4,ind3); %���,��ȥ���Լ����Լ�����߶�
    
    if  length(ind5)==2  %���������ٱ�,
        %--��ն������������
        closez=[closez startz];
        closex=[closex startx];
        %----�޸ĵ����Զ���η�Χ
        ind6=(dotarray(1,:)==startz & dotarray(2,:)==startx);
        POLYFLAG(ind6)=1;
        
        tempzx=linearray(:,ind5); %���������ٱߵ�������������
        %-----------������������������
        if ismember(ind5(1),ind1) %��ind5�е�һ��Ԫ�� �� ind1��
            tz1=tempzx(3,1); %ind1->(3,4)
            tx1=tempzx(4,1);%ind1
        else  %��ind5�е�һ��Ԫ�� �� ind2��
            tz1=tempzx(1,1);%ind2->(1,2)
            tx1=tempzx(2,1);
        end
        if ismember(ind5(2),ind1) %��ind5�еڶ���Ԫ�� �� ind1��
            tz2=tempzx(3,2);%ind1->(3,4)
            tx2=tempzx(4,2);
        else  %��ind5�е�һ��Ԫ�� �� ind2��
            tz2=tempzx(1,2);%ind2->(1,2)
            tx2=tempzx(2,2);
        end
        %------�պ�·�߷���ѡ��
        if tz1==tmz && tx1==tmx %���ĸ�����֮ǰ�ĵ㣬
            nextz=tz2;
            nextx=tx2;
        elseif tz2==tmz && tx2==tmx
            nextz=tz1;
            nextx=tx1;
        end
        %----ǰһ����
        tmz=startz;
        tmx=startx;
        %---��һ����
        startz=nextz;
        startx=nextx;
        
    else
        ind6=(dotarray(1,:)==startz & dotarray(2,:)==startx);
        POLYFLAG(ind6)=2; %�����������ٱߵģ�Ҳ�ų�����
        return;%ֱ�ӽ�������
    end % �Ƿ�������
end %while

if (startz==endz0 && startx==endx0) && length(ind5)==2 %���߷��
    %--��ն�������һ������
    closez=[closez endz0];
    closex=[closex endx0];
    %----�޸ĵ����Զ���η�Χ
    ind6=(dotarray(1,:)==endz0 & dotarray(2,:)==endx0);
    POLYFLAG(ind6)=1;
    
    closeflag=1;
end

end