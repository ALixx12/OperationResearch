% 该程序是求解产销平衡问题初始基可行解的西北角法
% 输入项
% Demand 为需求量或销量(列向量)；Supply为供应量或产地产量(列向量)
% 输出项
% X 为运输分配方案;b 记录了数字格(1 表示对应的是基变量，数字格，0 表示非基变量)
function [X,b]=NorthWest(Supply,Demand)
m = length(Supply);n=length(Demand);
i=1;j=1;
X=zeros(m,n);b=zeros(m,n);
while((i<=m)&&(j<=n))
    if Supply(i)<Demand(j)
        X(i,j)=Supply(i);
        b(i,j)=1;
        Demand(j)=Demand(j)-Supply(i);
        i=i+1;
    else
        X(i,j)=Demand(j);
        b(i,j)=1;
        Supply(i)=Supply(i)-Demand(j);
        j=j+1;
    end
end
