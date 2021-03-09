%limit equation 1
%验证不出来 差分Fourier谱
% clc;
% clear;
%format long;
function uw=limit(T)

%grids
J=128*2*16;
N=100000;
lada=1;%the coefficient of the nonlinear term
a=-128;b=128;

h=(b-a)/J;
T=T;
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
 u0(k,1)=exp(-(a+(k-1)*h)^2)/sqrt(pi);
        v0(k,1)=1/2*sech((a+(k-1)*h)^2)*sin(a+(k-1)*h);

    end
 %limit equation's initial value 
    u(:,1)=u0-i*v0;
    v(:,1)=conj(u0)-i*conj(v0);
%     FSU(:,1)=1/2*(u(:,1)+conj(v(:,1)));
%matrix
for k=1:J
    if k>J/2
    B(k)=-(k-J-1)^2;
    else
    B(k)=-(k-1)^2;
    end
end
B=B*(2*pi/(b-a))^2;

for n=1:N
ee=1;
%iteration value
u1=u(:,1);%take the last level as iterative initial value
v1=v(:,1);
u2=fft(u(:,1));
v2=fft(v(:,1));
while ee > 10^(-8)
    %nonlinear term
    unf1=fft(lada/32*(abs(u(:,1)).^2+abs(u1).^2+2*abs(v(:,1)).^2+2*abs(v1).^2).*(u(:,1)+u1));
    vnf1=fft(lada/32*(2*abs(u(:,1)).^2+2*abs(u1).^2+abs(v(:,1)).^2+abs(v1).^2).*(v(:,1)+v1));
    for l=1:J
        u(l,2)=(unf1(l)+(i/tau+B(l)/4)*u2(l))/(i/tau-B(l)/4);
        v(l,2)=(vnf1(l)+(i/tau+B(l)/4)*v2(l))/(i/tau-B(l)/4);
    end
    u(:,2)=ifft(u(:,2));
    v(:,2)=ifft(v(:,2));
    ee=sqrt(h)*norm(u(:,2)-u1,'fro');
    u1=u(:,2);
    v1=v(:,2);
end%while
u(:,1)=u(:,2);
v(:,1)=v(:,2);
if mod(n,1000)==0
    (n)/1000
    Us(:,(n)/1000)=u(:,2);
    tim((n)/1000)=(n)*tau;
    Vs(:,(n)/1000)=u(:,2);
end

end%time
for fe=7:7
    epsilon=10^(-fe);
    for fn=1:100
 FSU(:,fn,fe)=1/2*(exp(i/epsilon^2*tim(fn))*Us(:,fn)+exp(-i/epsilon^2*tim(fn))*conj(Vs(:,fn)));
    end
 %FSU1(:,fe)=1/2*(exp(i/(epsilon^2+1/2)*T)*u(:,2)+exp(-i/(epsilon^2+1/2)*T)*conj(v(:,2)));
 %x=linspace(a,b-h,J);
%  t=linspace(0,T,N+1);
%  [X,T]=meshgrid(t,x);
%  mesh(X,T,real(FSU))
%
% error(fe)=sqrt(h)*norm((FU(:,fe)-FSU(:,fe))/epsilon^0,'fro');
% errorf(fe)=norm((FU(:,fe)-FSU(:,fe))/epsilon^0,'inf');
end%epsilon
uw=FSU;
end
%对比


