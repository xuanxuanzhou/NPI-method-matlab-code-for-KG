%limit equation 3
%分裂步fourier方法求解极限方程
clc;clear;
%format long;
tic
%grids
J=2^8;
N=100000;
lada=1;%the coefficient of the nonlinear term
a=-128;b=128;

h=(b-a)/J;
T=10;
tau=T/N;
%uu=zeros(J,1000);

%the first level value of the numerical solution
u=zeros(J,2);
v=u;
u0=zeros(J,1);%一次性使用
v0=u0;
    for k=1:J
        %u0->\psi_0,v0->\psi_1
        %example one
%         u0(k)=1/2*(cos(3*(a+(k-1)*h))^2*sin(2*(a+(k-1)*h)))/(2-cos(a+(k-1)*h));
%         v0(k)=1/2*(cos(2*(a+(k-1)*h))*sin(a+(k-1)*h))/(2-cos(a+(k-1)*h));
        %example two 
%         u0(k)= (2+i)/sqrt(5)*cos((a+(k-1)*h));
%         v0(k)=(1+i)/sqrt(2)*sin((a+(k-1)*h))+1/2*cos((a+(k-1)*h));%一阶导数初值;
%         %example three
 u0(k,1)=exp(-(a+(k-1)*h)^2)/sqrt(pi);
        u1(k,1)=1/2*sech((a+(k-1)*h)^2)*sin(a+(k-1)*h);
    end
 %limit equation's initial value 
    u(:,1)=u0-i*v0;
    v(:,1)=conj(u0)-i*conj(v0);
    u0=abs(u(:,1)).^2;
    v0=abs(v(:,1)).^2;
    
%     FSU(:,1)=1/2*(u(:,1)+conj(v(:,1)));
%matrix
for k=1:J
    if k>J/2
    B(k)=-(k-J-1)^2;
    else
    B(k)=-(k-1)^2;
    end
end

for n=1:N


u2=fft(u(:,1));
v2=fft(v(:,1));
%第一步
    for l=1:J
        u(l,2)=(2*i/tau+B(l)/4)*u2(l)/(2*i/tau-B(l)/4);
        v(l,2)=(2*i/tau+B(l)/4)*v2(l)/(2*i/tau-B(l)/4);
    end
    u(:,2)=ifft(u(:,2));
    v(:,2)=ifft(v(:,2));
    %第二步
    for l=1:J
        u(l,2)=exp(-i*tau*1/8*(u0(l)+2*v0(l)))*u(l,2);
        v(l,2)=exp(-i*tau*1/8*(v0(l)+2*u0(l)))*v(l,2);
    end
    %第三步
 u2=fft(u(:,2));
v2=fft(v(:,2));
    for l=1:J
        u(l,2)=(2*i/tau+B(l)/4)*u2(l)/(2*i/tau-B(l)/4);
        v(l,2)=(2*i/tau+B(l)/4)*v2(l)/(2*i/tau-B(l)/4);
    end
    u(:,2)=ifft(u(:,2));
    v(:,2)=ifft(v(:,2));
    u(:,1)=u(:,2);
    v(:,1)=v(:,2);
%  FSU(:,n+1)=1/2*(exp(i/epsilon^2*n*tau)*u(:,2)+exp(-i/epsilon^2*n*tau)*conj(v(:,2)));
end%time
for fe=6
    epsilon=10^(-(fe-1));
 FSU(:,fe)=1/2*(exp(i/epsilon^2*T)*u(:,2)+exp(-i/epsilon^2*T)*conj(v(:,2)));
 %FSU1(:,fe)=1/2*(exp(i/(epsilon^2+1/2)*T)*u(:,2)+exp(-i/(epsilon^2+1/2)*T)*conj(v(:,2)));
% x=linspace(a,b-h,J);
% t=linspace(0,T,N+1);
% [X,T]=meshgrid(t,x);
% mesh(X,T,real(FSU))
% figure
% plot(x,real(FSU(:,N+1)-uN1))
% error(fe)=sqrt(h)*norm((FU(:,fe)-FSU(:,fe))/epsilon^0,'fro');
% errorf(fe)=norm((FU(:,fe)-FSU(:,fe))/epsilon^0,'inf');
end%epsilon
toc
%对比


