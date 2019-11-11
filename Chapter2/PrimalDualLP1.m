% ԭʼ��ż�ڵ��㷨������Թ滮
% min?c^T'x
% s.t.Ax=b
% Ax=b
% x>=0
% ���ż����Ϊ
% max b^Ty
% z=c�\A^Ty
% z>=0
% ��������Ҫ�����ʼ���е�
% ���� gamma ӦС�� 1 �������ܽӽ� 1����ȡ 0.9999 �ȣ�varepsilon Ϊ����Ҫ��
% mu ��Ӧ���䲻�ϼ�С��Ҳ��ȡ x^Tz/n
% ������� xstar,fstar �ֱ���ԭ��������Ž������Ŀ�꺯��ֵ
% (ystar,zstar),omegastar �ֱ��Ƕ�ż��������Ž������Ŀ�꺯��ֵ
function[xstar,ystar,zstar,fstar,omegastar]=PrimalDualLP1(A0,b0,c0,gamma,varepsilon)
tic;
[m0,n0]=size(A0);
x0=ones(n0,1);y0=zeros(m0,1);z0=ones(n0,1);
% ����� xi �� eta Ϊ�����˹�ԭʼ������˹���ż�������������������
% xi>(b�\A*x0)'*y0;eta>(A'*y0+z0�\c)'*x0.
theta=10;k=0;
xi=theta+(b0-A0*x0)'*y0;eta=theta+(A0'*y0+z0-c0)'*x0;
xstar=[x0;1;eta-(A0'*y0+z0-c0)'*x0];ystar=[y0;-1];zstar=[z0;xi-(b0-A0*x0)'*y0;1];
A=[A0,b0-A0*x0,zeros(m0,1);(A0'*y0+z0-c0)',0,1];b=[b0;eta];c=[c0;xi;0];
[m,n]=size(A);
while (abs(xstar(n-1))>=varepsilon)|(abs(ystar(m))>=varepsilon)
    xi=theta+(b0-A0*x0)'*y0;eta=theta+(A0'*y0+z0-c0)'*x0;
    x=[x0;1;eta-(A0'*y0+z0-c0)'*x0];y=[y0;-1];z=[z0;xi-(b0-A0*x0)'*y0;1];
    A=[A0,b0-A0*x0,zeros(m0,1);(A0'*y0+z0-c0)',0,1];b=[b0;eta];c=[c0;xi;0];
    while x'*z>varepsilon
        Z=diag(z);X=diag(x);beta=1-0.01/sqrt(n);mu=x'*z/n;MuPlus=beta*mu;
        A1=[Z,zeros(n,m),X;A,zeros(m,m),zeros(m,n);zeros(n,n),A',eye(n,n)];
        b1=[(MuPlus*eye(n,n)-X*Z)*ones(n,1);zeros(m,1);zeros(n,1)];
        %����Ϊ�ⷽ�����ѡ����Ԫ Gauss ��Ԫ����ͨ���ⷽ����õ� x,y,z ��������
        n1=length(b1);
        A2=[A1,b1];
        for k=1:n1-1 % ѡ��Ԫ
            [Ap,p]=max(abs(A2(k:n1,k)));
            p=p+k-1;
            if p>k
                t=A2(k,:);A2(k,:)=A2(p,:);A2(p,:)=t;
            end
            %����Ϊ��Ԫ
            A2((k+1):n1,(k+1):(n1+1))=A2((k+1):n1,(k+1):(n1+1))-A2((k+1):n1,k)/A2(k,k)*A2(k,(k+ 1):(n1+1));
            A2((k+1):n1,k)=zeros(n1-k,1);
        end
        % ����Ϊ�ش�
        x1=zeros(n1,1);
        x1(n1)=A2(n1,n1+1)/A2(n1,n1);
        for k=n1-1:-1:1
            x1(k,:)=(A2(k,n1+1)-A2(k,(k+1):n1)*x1((k+1):n1))/A2(k,k);
        end
        Solve=x1;
        %���ϳ����Ҳ��ͨ��Solve=Gauss2(A1,b1,1);������ѡ����ԪGauss��Ԫ���ⷽ����
        dx=Solve(1:n);dy=Solve(n+1:n+m);dz=Solve(n+m+1:2*n+m);
        k=0;j=0;
        for i=1:n
            if dx(i)<0
                k=k+1;
                AlphaP(k)=-(x(i)/dx(i));
            end
            if dz(i)<0
                j=j+1;
                AlphaD(j)=-(z(i)/dz(i));
            end
        end
        AlphaP=min(AlphaP);AlphaD=min(AlphaD);
        AlphaMax=min(AlphaP,AlphaD);
        alpha=gamma*AlphaMax;
        x=x+alpha*dx;y=y+alpha*dy;z=z+alpha*dz;
    end
    xstar=x;ystar=y;zstar=z;
    %if (abs(xstar(n-1))<varepsilon)&(abs(ystar(m))<varepsilon)
    %   x=xstar(1:n�\2);y=ystar(1:m�\1);z=zstar(1:n�\2);
    %end
    k=k+1;theta=10*theta;
end
x=xstar(1:n-2);y=ystar(1:m-1);z=zstar(1:n-2);
disp('ԭ�������Ž�Ϊ xstar,��ż��������Ž� ystar,Ŀ�꺯���ֱ�Ϊ fstar,omegastar��')
xstar=x;ystar=y;zstar=z;fstar=c0'*xstar;omegastar=b0'*ystar;k,theta
toc;