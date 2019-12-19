%单纯形法求解目标规划问题，采用分层、逐次利用单纯形法求解
%须将目标规划化为标准形式
%输入项：A 为约束方程组系数矩阵；b 为约束方程右端项
%c 为达成函数依优先目标按行输入，即c 的第一行是第一优先目标不考虑优先因子时的目标函数
%第二行是只考虑第二优先目标(不考虑优先因子)时的目标函数
%以下类推
function x=GoalProg(A,b,c)
[m0,n0]=size(A);[m1,n1]=size(c);
if n1~=n0
    disp('目标函数矩阵输入有误')
    return
end
for i=1:m1
    [xstar,fxstar,A0,IB,iter]=MMSimplex(A,b,-c(i,:)');
    A(m0+i,:)=c(i,:);b(m0+i)=-fxstar;
end
x=xstar;