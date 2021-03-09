%limit equation  for quadratic nonlinear term
%验证不出来 差分Fourier谱
% clc;
% clear;
%format long;
function uw=limitq(T)

%grids
J=32*2*16;
N=100000;
lada=1;%the coefficient of the nonlinear term
a=-32;b=32;

h=(b-a)/J;
T=T;
tau=T/N;
%uu=zeros(J,1000);

%the first level value of the numerical solution
u=zeros(J,2);
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
    u2=fft(u(:,1));
    for l=1:J
        u(l,2)=(i/tau+B(l)/4)*u2(l)/(i/tau-B(l)/4);
    end
    u(:,2)=ifft(u(:,2));
    u(:,1)=u(:,2);
    
    if mod(n,1000)==0
        (n)/1000
        Us(:,(n)/1000)=u(:,2);
        tim((n)/1000)=(n)*tau;
    end
    
end%time
for fe=1:4
    epsilon=10^(-fe);
    for fn=1:100
        FSU(:,fn,fe)=1/2*(exp(i/epsilon^2*tim(fn))*Us(:,fn)+exp(-i/epsilon^2*tim(fn))*conj(Us(:,fn)));
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


