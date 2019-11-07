% Ssimplex.m  利用单纯形法的 MATLAB 程序
% max c'x 
% a.t. 
% ax=b 
% x>=0 
%  这里 A\in R^{m\times n},c,x'\im R^m,b\geq 0 
%  且矩阵 A 中有一个单位子矩阵，不需要引入人工变量
function [xstar,fxstar,iter]=Ssimplex(A,b,c) 
[m,n]=size(A);E=eye(m);IB=zeros(1,m);k=0; 
for i = 1:m 
        for j=1:n 
                if A(:,j)==E(:,i) 
                        IB(i)=j;SA(i)=j;%IB 记录基变量下标，SA 记录松弛变量下标
                elseif A(:,j)==(-E(:,i))
                        SA(i)=j;%SA 也记录剩余变量（松弛变量）下标
                end 
        end 
end 
A0=[b,A];N=1:n;N(IB)=[];IN=N;x(IB)=A0(:,1)'; 
x(IN)=zeros(1,length(IN));cB=c(IB); 
%IN 为非基变量下标
sigma=c'-cB'*A0(:,2:n+1);t=length(find(sigma>0)); 
%计算原问题的检验数并假设检验数中有 t 个大于零的检验数
while t~=0 
        [sigmaJ,jj]=max(sigma); 
    %这里的 jj 是 sigma 中值最大者所在列，即 A0 中的第 jj+1 列（A0 中的第一列为 b），该列对应的非基变量 x(jj)为换入变量，而 sigmaJ 则是相应的检验数
        tt=find(A0(:,jj+1)>0);kk=length(tt); 
    %检查增广系数矩阵 A0 中第 jj+1 列元素是否有大于零的元素
        if kk==0 
                disp('原问题为无界解') 
        else 
                theta =zeros(1,kk); 
                for i=1:kk 
                        theta(i)=A0(tt(i),1)/A0(tt(i),jj+1); 
                end 
                [thetaI,ii]=min(theta);Temp=tt(ii); 
%比值最小的 theta 值， 选择换出变量。 这时 A0(Temp,jj+1)为旋转主元
                for i=1:m 
                        if i~=Temp 
                                A0(i,:)=A0(i,:)-(A0(Temp,:)/A0(Temp,jj+1))*A0(i,jj+1); 
                        else 
                                A0(Temp,:)=A0(Temp,:)/A0(Temp,jj+1); 
                        end 
                end 
                TT=IB(Temp);IB(Temp)=jj;IN(jj)=TT;x(IB)=A(:,1)'; 
                N=1:n;N(IB)=[];IN=N;x(IN)=zeros(1,length(IN));cB=c(IB);
        %新的基可行解及其价值系数
                sigma=c'-cB'*A0(:,2:n+1);t=length(find(sigma>0)); 
        %再次计算检验数并假设检验数中有 t 个大于零的检验数
        end 
        k=k+1; 
end 
IB 
IN 
B=A(:,IB); 
%InverseOfB=inv(B) 
%这是基矩阵 B 的逆矩阵，用于灵敏度分析。若不做灵敏度分析，则将其注释掉
xstar=x;fxstar=x(IB)*c(IB);iter=k; 