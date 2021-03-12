%% program handbook
%没有给算子D，exp(aD)近似
%the program for the Klein-Gordon equation in the non-relativistic limit
%regime
%sinF.m is same as sinFu.m(replace the D as sin (ttD)/tt)
%Nested Picard method with finite difference scheme/ Fourier pseudospectral method
%% integral end start computing
%parameters
% tic;
% clc;clear;
function uw=sinFu(T)

T=T;a=-128;b=128;
lada=1;%the coefficient of nonlinear term
%符号计算
syms s w epsilo real
f=exp(1i*(w-s)/epsilo^2)-exp(-1i*(w-s)/epsilo^2);%cf=-f
g=exp(1i*(w-s)/epsilo^2)+exp(-1i*(w-s)/epsilo^2);%cg=g
%first-order used funciton
f1{1}=symfun(int(f*exp(3*1i*s/epsilo^2),s,0,w),[s,epsilo]);
f1{2}=symfun(int(f*exp(1i*s/epsilo^2),s,0,w),[s,epsilo]);
f1{3}=symfun(int(f*exp(-1i*s/epsilo^2),s,0,w),[s,epsilo]);
f1{4}=symfun(int(f*exp(-3*1i*s/epsilo^2),s,0,w),[s,epsilo]);
f1{1}=subs(f1{1},w,s);
f1{2}=subs(f1{2},w,s);
f1{3}=subs(f1{3},w,s);
f1{4}=subs(f1{4},w,s);
g1{1}=symfun(int(g*exp(3*1i*s/epsilo^2),s,0,w),[s,epsilo]);
g1{2}=symfun(int(g*exp(1i*s/epsilo^2),s,0,w),[s,epsilo]);
g1{3}=symfun(int(g*exp(-1i*s/epsilo^2),s,0,w),[s,epsilo]);
g1{4}=symfun(int(g*exp(-3*1i*s/epsilo^2),s,0,w),[s,epsilo]);
g1{1}=subs(g1{1},w,s);
g1{2}=subs(g1{2},w,s);
g1{3}=subs(g1{3},w,s);
g1{4}=subs(g1{4},w,s);

ff1{1}=symfun(int(s*f*exp(3*1i*s/epsilo^2),s,0,w),[s,epsilo]);
ff1{2}=symfun(int(s*f*exp(1i*s/epsilo^2),s,0,w),[s,epsilo]);
ff1{3}=symfun(int(s*f*exp(-1i*s/epsilo^2),s,0,w),[s,epsilo]);
ff1{4}=symfun(int(s*f*exp(-3*1i*s/epsilo^2),s,0,w),[s,epsilo]);
ff1{1}=subs(ff1{1},w,s);
ff1{2}=subs(ff1{2},w,s);
ff1{3}=subs(ff1{3},w,s);
ff1{4}=subs(ff1{4},w,s);
gg1{1}=symfun(int(s*g*exp(3*1i*s/epsilo^2),s,0,w),[s,epsilo]);
gg1{2}=symfun(int(s*g*exp(1i*s/epsilo^2),s,0,w),[s,epsilo]);
gg1{3}=symfun(int(s*g*exp(-1i*s/epsilo^2),s,0,w),[s,epsilo]);
gg1{4}=symfun(int(s*g*exp(-3*1i*s/epsilo^2),s,0,w),[s,epsilo]);
gg1{1}=subs(gg1{1},w,s);
gg1{2}=subs(gg1{2},w,s);
gg1{3}=subs(gg1{3},w,s);
gg1{4}=subs(gg1{4},w,s);

fs2{1}=symfun(int(s*s*f*exp(3*1i*s/epsilo^2),s,0,w),[s,epsilo]);
fs2{2}=symfun(int(s*s*f*exp(1i*s/epsilo^2),s,0,w),[s,epsilo]);
fs2{3}=symfun(int(s*s*f*exp(-1i*s/epsilo^2),s,0,w),[s,epsilo]);
fs2{4}=symfun(int(s*s*f*exp(-3*1i*s/epsilo^2),s,0,w),[s,epsilo]);
fs2{1}=subs(fs2{1},w,s);
fs2{2}=subs(fs2{2},w,s);
fs2{3}=subs(fs2{3},w,s);
fs2{4}=subs(fs2{4},w,s);
gs2{1}=symfun(int(s*s*g*exp(3*1i*s/epsilo^2),s,0,w),[s,epsilo]);
gs2{2}=symfun(int(s*s*g*exp(1i*s/epsilo^2),s,0,w),[s,epsilo]);
gs2{3}=symfun(int(s*s*g*exp(-1i*s/epsilo^2),s,0,w),[s,epsilo]);
gs2{4}=symfun(int(s*s*g*exp(-3*1i*s/epsilo^2),s,0,w),[s,epsilo]);
gs2{1}=subs(gs2{1},w,s);
gs2{2}=subs(gs2{2},w,s);
gs2{3}=subs(gs2{3},w,s);
gs2{4}=subs(gs2{4},w,s);
%right
% %second-order used function
sf1{1}=symfun(int((w-s)*g*exp(3*1i*s/epsilo^2),s,0,w),[s,epsilo]);
sf1{2}=symfun(int((w-s)*g*exp(1i*s/epsilo^2),s,0,w),[s,epsilo]);
sf1{3}=symfun(int((w-s)*g*exp(-1i*s/epsilo^2),s,0,w),[s,epsilo]);
sf1{4}=symfun(int((w-s)*g*exp(-3*1i*s/epsilo^2),s,0,w),[s,epsilo]);
sf1{1}=subs(sf1{1},w,s);
sf1{2}=subs(sf1{2},w,s);
sf1{3}=subs(sf1{3},w,s);
sf1{4}=subs(sf1{4},w,s);
sg1{1}=symfun(int((w-s)*f*exp(3*1i*s/epsilo^2),s,0,w),[s,epsilo]);
sg1{2}=symfun(int((w-s)*f*exp(1i*s/epsilo^2),s,0,w),[s,epsilo]);
sg1{3}=symfun(int((w-s)*f*exp(-1i*s/epsilo^2),s,0,w),[s,epsilo]);
sg1{4}=symfun(int((w-s)*f*exp(-3*1i*s/epsilo^2),s,0,w),[s,epsilo]);
sg1{1}=subs(sg1{1},w,s);
sg1{2}=subs(sg1{2},w,s);
sg1{3}=subs(sg1{3},w,s);
sg1{4}=subs(sg1{4},w,s);
% 

fff1{1}=symfun(int((w-s)*s*g*exp(3*1i*s/epsilo^2),s,0,w),[s,epsilo]);
fff1{2}=symfun(int((w-s)*s*g*exp(1i*s/epsilo^2),s,0,w),[s,epsilo]);
fff1{3}=symfun(int((w-s)*s*g*exp(-1i*s/epsilo^2),s,0,w),[s,epsilo]);
fff1{4}=symfun(int((w-s)*s*g*exp(-3*1i*s/epsilo^2),s,0,w),[s,epsilo]);
fff1{1}=subs(fff1{1},w,s);
fff1{2}=subs(fff1{2},w,s);
fff1{3}=subs(fff1{3},w,s);
fff1{4}=subs(fff1{4},w,s);
ggg1{1}=symfun(int((w-s)*s*f*exp(3*1i*s/epsilo^2),s,0,w),[s,epsilo]);
ggg1{2}=symfun(int((w-s)*s*f*exp(1i*s/epsilo^2),s,0,w),[s,epsilo]);
ggg1{3}=symfun(int((w-s)*s*f*exp(-1i*s/epsilo^2),s,0,w),[s,epsilo]);
ggg1{4}=symfun(int((w-s)*s*f*exp(-3*1i*s/epsilo^2),s,0,w),[s,epsilo]);
ggg1{1}=subs(ggg1{1},w,s);
ggg1{2}=subs(ggg1{2},w,s);
ggg1{3}=subs(ggg1{3},w,s);
ggg1{4}=subs(ggg1{4},w,s);

% 
fws2{1}=symfun(int((w-s)*(w-s)*f*exp(3*1i*s/epsilo^2),s,0,w),[s,epsilo]);
fws2{2}=symfun(int((w-s)*(w-s)*f*exp(1i*s/epsilo^2),s,0,w),[s,epsilo]);
fws2{3}=symfun(int((w-s)*(w-s)*f*exp(-1i*s/epsilo^2),s,0,w),[s,epsilo]);
fws2{4}=symfun(int((w-s)*(w-s)*f*exp(-3*1i*s/epsilo^2),s,0,w),[s,epsilo]);
fws2{1}=subs(fws2{1},w,s);
fws2{2}=subs(fws2{2},w,s);
fws2{3}=subs(fws2{3},w,s);
fws2{4}=subs(fws2{4},w,s);
gws2{1}=symfun(int((w-s)*(w-s)*g*exp(3*1i*s/epsilo^2),s,0,w),[s,epsilo]);
gws2{2}=symfun(int((w-s)*(w-s)*g*exp(1i*s/epsilo^2),s,0,w),[s,epsilo]);
gws2{3}=symfun(int((w-s)*(w-s)*g*exp(-1i*s/epsilo^2),s,0,w),[s,epsilo]);
gws2{4}=symfun(int((w-s)*(w-s)*g*exp(-3*1i*s/epsilo^2),s,0,w),[s,epsilo]);
gws2{1}=subs(gws2{1},w,s);
gws2{2}=subs(gws2{2},w,s);
gws2{3}=subs(gws2{3},w,s);
gws2{4}=subs(gws2{4},w,s);
% 

% secondorder n.2(2)
cf1{1}=symfun(int(-f*exp(-3*1i*s/epsilo^2),s,0,w),[s,epsilo]);
cf1{2}=symfun(int(-f*exp(-1i*s/epsilo^2),s,0,w),[s,epsilo]);
cf1{3}=symfun(int(-f*exp(1i*s/epsilo^2),s,0,w),[s,epsilo]);
cf1{4}=symfun(int(-f*exp(3*1i*s/epsilo^2),s,0,w),[s,epsilo]);
cf1{1}=subs(cf1{1},w,s);
cf1{2}=subs(cf1{2},w,s);
cf1{3}=subs(cf1{3},w,s);
cf1{4}=subs(cf1{4},w,s);

f2{1}=symfun(exp(2*1i*s/epsilo^2),[s,epsilo]);
f2{2}=symfun(exp(0*1i*s/epsilo^2),[s,epsilo]);
f2{3}=symfun(exp(-2*1i*s/epsilo^2),[s,epsilo]);
%二阶拓展
for k=1:3
    for kk=1:4
        f5{k,kk}=symfun(int(f*f2{k}*f1{kk},s,0,w),[s,epsilo]);
        f5{k,kk}=subs(f5{k,kk},w,s);
        f6{k,kk}=symfun(int(f*f2{k}*cf1{kk},s,0,w),[s,epsilo]);
        f6{k,kk}=subs(f6{k,kk},w,s);
        fs5{k,kk}=symfun(int(s*f*f2{k}*f1{kk},s,0,w),[s,epsilo]);
        fs5{k,kk}=subs(fs5{k,kk},w,s);
        fs6{k,kk}=symfun(int(s*f*f2{k}*cf1{kk},s,0,w),[s,epsilo]);
        fs6{k,kk}=subs(fs6{k,kk},w,s);

        f7{k,kk}=symfun(int(g*f2{k}*f1{kk},s,0,w),[s,epsilo]);
        f7{k,kk}=subs(f7{k,kk},w,s);
        f8{k,kk}=symfun(int(g*f2{k}*cf1{kk},s,0,w),[s,epsilo]);
        f8{k,kk}=subs(f8{k,kk},w,s);
         fs7{k,kk}=symfun(int(s*g*f2{k}*f1{kk},s,0,w),[s,epsilo]);
        fs7{k,kk}=subs(fs7{k,kk},w,s);
        fs8{k,kk}=symfun(int(s*g*f2{k}*cf1{kk},s,0,w),[s,epsilo]);
        fs8{k,kk}=subs(fs8{k,kk},w,s);

        fws5{k,kk}=symfun(int((w-s)*g*f2{k}*f1{kk},s,0,w),[s,epsilo]);
        fws5{k,kk}=subs(fws5{k,kk},w,s);
        fws6{k,kk}=symfun(int((w-s)*g*f2{k}*cf1{kk},s,0,w),[s,epsilo]);
        fws6{k,kk}=subs(fws6{k,kk},w,s);
        fws7{k,kk}=symfun(int((w-s)*f*f2{k}*f1{kk},s,0,w),[s,epsilo]);
        fws7{k,kk}=subs(fws7{k,kk},w,s);
        fws8{k,kk}=symfun(int((w-s)*f*f2{k}*cf1{kk},s,0,w),[s,epsilo]);
        fws8{k,kk}=subs(fws8{k,kk},w,s);
    end
end
%third order
%others
f9{1}=symfun(exp(1i*s/epsilo^2),[s,epsilo]);
f9{2}=symfun(exp(-1i*s/epsilo^2),[s,epsilo]);


for k=1:4
    for kk=1:4
    for kkk=1:2
        f10{k,kk,kkk}=symfun(int(f*f1{k}*cf1{kk}*f9{kkk},s,0,w),[s,epsilo]);
        f10{k,kk,kkk}=subs(f10{k,kk,kkk},w,s);

        f12{k,kk,kkk}=symfun(int(f*f1{k}*f1{kk}*f9{kkk},s,0,w),[s,epsilo]);
        f12{k,kk,kkk}=subs(f12{k,kk,kkk},w,s);

        g10{k,kk,kkk}=symfun(int(g*f1{k}*cf1{kk}*f9{kkk},s,0,w),[s,epsilo]);
        g10{k,kk,kkk}=subs(g10{k,kk,kkk},w,s);

        g12{k,kk,kkk}=symfun(int(g*f1{k}*f1{kk}*f9{kkk},s,0,w),[s,epsilo]);
        g12{k,kk,kkk}=subs(g12{k,kk,kkk},w,s);
    end
    end
end

% %第三项
for kkk=1:3
for k=1:3
    for kk=1:4
        %f
        f16{kkk,k,kk}=symfun(int(f*f2{kkk}*f5{k,kk},s,0,w),[s,epsilo]);
        f16{kkk,k,kk}=subs(f16{kkk,k,kk},w,s);
        f17{kkk,k,kk}=symfun(int(f*f2{kkk}*conj(f5{k,kk}),s,0,w),[s,epsilo]);
        f17{kkk,k,kk}=subs(f17{kkk,k,kk},w,s);
        f18{kkk,k,kk}=symfun(int(f*f2{kkk}*f6{k,kk},s,0,w),[s,epsilo]);
        f18{kkk,k,kk}=subs(f18{kkk,k,kk},w,s);
        f19{kkk,k,kk}=symfun(int(f*f2{kkk}*conj(f6{k,kk}),s,0,w),[s,epsilo]);
        f19{kkk,k,kk}=subs(f19{kkk,k,kk},w,s);
        f20{kkk,k,kk}=symfun(int(g*f2{kkk}*f5{k,kk},s,0,w),[s,epsilo]);
        f20{kkk,k,kk}=subs(f20{kkk,k,kk},w,s);
        f21{kkk,k,kk}=symfun(int(g*f2{kkk}*conj(f5{k,kk}),s,0,w),[s,epsilo]);
        f21{kkk,k,kk}=subs(f21{kkk,k,kk},w,s);
        f22{kkk,k,kk}=symfun(int(g*f2{kkk}*f6{k,kk},s,0,w),[s,epsilo]);
        f22{kkk,k,kk}=subs(f22{kkk,k,kk},w,s);
        f23{kkk,k,kk}=symfun(int(g*f2{kkk}*conj(f6{k,kk}),s,0,w),[s,epsilo]);
        f23{kkk,k,kk}=subs(f23{kkk,k,kk},w,s);
    end
end
end
for k=1:3
    for kk=1:4
        f24{k,kk}=symfun(int(f*f2{k}*ff1{kk},s,0,w),[s,epsilo]);
        f24{k,kk}=subs(f24{k,kk},w,s);
        f25{k,kk}=symfun(int(f*f2{k}*conj(ff1{kk}),s,0,w),[s,epsilo]);
        f25{k,kk}=subs(f25{k,kk},w,s);
        f26{k,kk}=symfun(int(g*f2{k}*ff1{kk},s,0,w),[s,epsilo]);
        f26{k,kk}=subs(f26{k,kk},w,s);
        f27{k,kk}=symfun(int(g*f2{k}*conj(ff1{kk}),s,0,w),[s,epsilo]);
        f27{k,kk}=subs(f27{k,kk},w,s);
        f28{k,kk}=symfun(int(f*f2{k}*sf1{kk},s,0,w),[s,epsilo]);
        f28{k,kk}=subs(f28{k,kk},w,s);
        f29{k,kk}=symfun(int(f*f2{k}*conj(sf1{kk}),s,0,w),[s,epsilo]);
        f29{k,kk}=subs(f29{k,kk},w,s);
        f30{k,kk}=symfun(int(g*f2{k}*sf1{kk},s,0,w),[s,epsilo]);
        f30{k,kk}=subs(f30{k,kk},w,s);
        f31{k,kk}=symfun(int(g*f2{k}*conj(sf1{kk}),s,0,w),[s,epsilo]);
        f31{k,kk}=subs(f31{k,kk},w,s);
    end
end
%% starting

for fe=1:4
    disp(fe);
    epsilon=10^(-fe);
for  fN=1:1
    disp(fN);
    
%M=2^3*2^(fN-1);
M=128*2*16;

h=(b-a)/M;
N=10000;
%u=zeros(M,N);
%N=20*2^(fN-1);
tau=T/N;
tt=tau;
for k=1:4
            f1a(k)=double(f1{k}(tau,epsilon));
        sf1a(k)=double(sf1{k}(tau,epsilon));
        ff1a(k)=double(ff1{k}(tau,epsilon));
        fws2a(k)=double(fws2{k}(tau,epsilon));
        fff1a(k)=double(fff1{k}(tau,epsilon));
        fs2a(k)=double(fs2{k}(tau,epsilon));
        g1a(k)=double(g1{k}(tau,epsilon));
        sg1a(k)=double(sg1{k}(tau,epsilon));
        gg1a(k)=double(gg1{k}(tau,epsilon));
        gws2a(k)=double(gws2{k}(tau,epsilon));
       ggg1a(k)=double(ggg1{k}(tau,epsilon));
        gs2a(k)=double(gs2{k}(tau,epsilon));
end
for k=1:3
    for kk=1:4
                  f5a(k,kk)=double(f5{k,kk}(tau,epsilon));
           f6a(k,kk)=double(f6{k,kk}(tau,epsilon));
           fs5a(k,kk)=double(fs5{k,kk}(tau,epsilon));
           fs6a(k,kk)=double(fs6{k,kk}(tau,epsilon));
           fws5a(k,kk)=double(fws5{k,kk}(tau,epsilon));
           fws6a(k,kk)=double(fws6{k,kk}(tau,epsilon));
           f7a(k,kk)=double(f7{k,kk}(tau,epsilon));
           f8a(k,kk)=double(f8{k,kk}(tau,epsilon));
           fs7a(k,kk)=double(fs7{k,kk}(tau,epsilon));
           fs8a(k,kk)=double(fs8{k,kk}(tau,epsilon));
           fws7a(k,kk)=double(fws7{k,kk}(tau,epsilon));
           fws8a(k,kk)=double(fws8{k,kk}(tau,epsilon));
    end
end
for k=1:4
    for kk=1:4
        for kkk=1:2
                           f10a(k,kk,kkk)=double(f10{k,kk,kkk}(tau,epsilon));
               f12a(k,kk,kkk)=double(f12{k,kk,kkk}(tau,epsilon));
                g10a(k,kk,kkk)=double(g10{k,kk,kkk}(tau,epsilon));
               g12a(k,kk,kkk)=double(g12{k,kk,kkk}(tau,epsilon));
        end
    end
end
for k=1:3
    for kk=1:3
        for kkk=1:4
             f16a(k,kk,kkk)=double(f16{k,kk,kkk}(tau,epsilon));
            f17a(k,kk,kkk)=double(f17{k,kk,kkk}(tau,epsilon));
            f18a(k,kk,kkk)=double(f18{k,kk,kkk}(tau,epsilon));
            f19a(k,kk,kkk)=double(f19{k,kk,kkk}(tau,epsilon));
            f20a(k,kk,kkk)=double(f20{k,kk,kkk}(tau,epsilon));
            f21a(k,kk,kkk)=double(f21{k,kk,kkk}(tau,epsilon));
            f22a(k,kk,kkk)=double(f22{k,kk,kkk}(tau,epsilon));
            f23a(k,kk,kkk)=double(f23{k,kk,kkk}(tau,epsilon));
        end
    end
end
for k=1:3
    for kk=1:4
         f24a(k,kk)=double(f24{k,kk}(tau,epsilon));
        f25a(k,kk)=double(f25{k,kk}(tau,epsilon));
        f26a(k,kk)=double(f26{k,kk}(tau,epsilon));
        f27a(k,kk)=double(f27{k,kk}(tau,epsilon));
        f28a(k,kk)=double(f28{k,kk}(tau,epsilon));
        f29a(k,kk)=double(f29{k,kk}(tau,epsilon));
        f30a(k,kk)=double(f30{k,kk}(tau,epsilon));
        f31a(k,kk)=double(f31{k,kk}(tau,epsilon));
    end
end
%initialization
u0=zeros(M,1);
u1=u0;
for k=1:M
    %the first set of initial conditions相似文章方程一致
%     UUU(k,1)=exp(i*7)*exp(i*(1+a+(k-1)*h));
%     u0(k,1)=exp(i*7)*exp(i*(a+(k-1)*h));
%     u1(k,1)=exp(i*7)*i*exp(i*(a+(k-1)*h));
%     u0(k,1)=1/2*(cos(3*(a+(k-1)*h))^2*sin(2*(a+(k-1)*h)))/(2-cos(a+(k-1)*h));
%     u1(k,1)=1/(2*epsilon^2)*(cos(2*(a+(k-1)*h))*sin(a+(k-1)*h))/(2-cos(a+(k-1)*h));%一阶导数初值
%     %the second set of initial conditions(variable
%     separation)klein-gordon方程有负号
%     u0(k,1)=(1+i)*exp(-(a+(k-1)*h)^2/2);
%     u1(k,1)=1/epsilon^2*3*exp(-(a+(k-1)*h)^2/2);%一阶导数初值
    %the third set of initial conditionsklein-gordon3方程一致lada=-1
%     u0(k,1)=(2+i)/sqrt(5)*cos(a+(k-1)*h);
%     u1(k,1)=1/epsilon^2*((1+i)/sqrt(2)*sin(a+(k-1)*h)+1/2*cos(a+(k-1)*h));%一阶导数初值
%     %the fourth set of initial conditions(the dirichlet
% %     boundary)baolaoshi方程有负号
%     u0(k,1)=exp(-(a+(k-1)*h)^2)/sqrt(pi);
%     u1(k,1)=1/epsilon^2*sech((a+(k-1)*h)^2)*sin(a+(k-1)*h);%一阶导数初值
%example three
%         u0(k,1)= exp(-4*(a+(k-1)*h)^2);
%         u1(k,1)=1/epsilon^2*u0(k);%一阶导数初值;
        %exact solution
%         u0(k,1)= exp(i*(a+(k-1)*h));%exp(i*(t/8+a+(k-1)*h))
%         u1(k,1)=1/epsilon^2*u0(k);%一阶导数初值;

 u0(k,1)=exp(-(a+(k-1)*h)^2)/sqrt(pi);
        u1(k,1)=1/epsilon^2*1/2*sech((a+(k-1)*h)^2)*sin(a+(k-1)*h);
end
% sptial matrix
mul=zeros(M,1);
betal=mul;
betag=mul;
for l=1:M
    %notice the solver of 0->0(the first method for l: 0~M/2-1<=>1~M/2, -M/2~-1<=>M/2+1~M)
        if l>M/2
            mul(l)=2*pi*(l-M-1)/(b-a);
        else
            mul(l)=2*pi*(l-1)/(b-a);
        end
        betal(l,1)=sqrt(1+epsilon^2*mul(l)^2)/(epsilon^2);
end
betag=betal-1/epsilon^2;
betag=sin(tt*(betag))/tt;
%% computing preparation
u0=fft(u0);
u1=fft(u1);
Y=lada./(2*1i*betal*epsilon^2);
z=lada/(2*epsilon^2);
%first-order auxillary variables
G1=zeros(M,4);G56=G1;G2=zeros(M,3);G4=zeros(M,3);
A=zeros(M,1);B=A;C=A;D=A;

%% computing
for tk=1:N %notice the tk can not be used again below
A=(1i*betal.*u0+u1)./(2*1i*betal);
B=(1i*betal.*u0-u1)./(2*1i*betal); 
C=(1i*betal.*u0+u1)/2; 
D= -1*(1i*betal.*u0-u1)/2; 
% zero order 
for l=1:M
u0(l)=exp(1i*tau*betal(l))*A(l)+exp(-1i*tau*betal(l))*B(l);
u1(l)=exp(1i*tau*betal(l))*C(l)+exp(-1i*tau*betal(l))*D(l);
end
E1=1i*betag.*A;
F1=-1i*betag.*B;
E1=ifft(E1);
F1=ifft(F1);
E2=-1/2*betag.^2.*A;
F2=-1/2*betag.^2.*B;
E2=ifft(E2);
F2=ifft(F2);
A=ifft(A);
B=ifft(B);

%G1物理空间
G1(:,1)=fft(A.^2.*conj(B));
G1(:,2)=fft(A.^2.*conj(A)+2*A.*B.*conj(B));
G1(:,3)=fft(conj(B).*B.^2+2*conj(A).*A.*B);
G1(:,4)=fft(conj(A).*B.^2);

G56(:,1)=fft(A.^2.*conj(F1))+2*fft(A.*conj(B).*E1);
G56(:,2)=fft(A.^2.*conj(E1)+2*A.*B.*conj(F1))+2*fft(A.*conj(B).*F1+(abs(A).^2+abs(B).^2).*E1);
G56(:,3)=fft(conj(F1).*B.^2+2*conj(E1).*A.*B)+2*fft(conj(A).*B.*E1+(abs(A).^2+abs(B).^2).*F1);
G56(:,4)=fft(conj(E1).*B.^2)+2*fft(conj(A).*B.*F1);

G11(:,1)=fft(A.^2.*conj(F2)+ E1.^2.*conj(B)+ 2*A.*E1.*conj(F1)+ 2*A.*conj(B).*E2);
G11(:,2)=fft(A.^2.*conj(E2)+2*A.*B.*conj(F2)+ E1.^2.*conj(A)+2*E1.*F1.*conj(B)+ 2*(B.*E1.*conj(F1)+(abs(E1).^2+abs(F1).^2).*A)+ 2*(A.*conj(B).*F2+(abs(A).^2+abs(B).^2).*E2));
G11(:,3)=fft(B.^2.*conj(F2)+2*A.*B.*conj(E2)+ F1.^2.*conj(B)+2*F1.*E1.*conj(A)+ 2*(A.*F1.*conj(E1)+(abs(E1).^2+abs(F1).^2).*B)+ 2*(B.*conj(A).*E2+(abs(A).^2+abs(B).^2).*F2));
G11(:,4)=fft(B.^2.*conj(E2)+ F1.^2.*conj(A)+ 2*B.*F1.*conj(E1)+ 2*B.*conj(A).*F2);
%first order
    %nonlinear term
   nontf1=zeros(M,1);
   nontg1=nontf1;
    for k=1:4%增加一行，阶数增加一阶
        
        nontf1=nontf1+f1a(k)*G1(:,k)...
         +sf1a(k)*i*betag.*G1(:,k)+ff1a(k)*G56(:,k)...
            +fws2a(k)*-1/2*betag.^2.*G1(:,k)+fff1a(k)*i*betag.*G56(:,k)+fs2a(k)*G11(:,k);
       
        nontg1=nontg1+g1a(k)*G1(:,k)...
         +sg1a(k)*i*betag.*G1(:,k)+gg1a(k)*G56(:,k)...
         +gws2a(k)*-1/2*betag.^2.*G1(:,k) +ggg1a(k)*i*betag.*G56(:,k)+gs2a(k)*G11(:,k);
    
    end

% second order
for k=1:4
    G1(:,k)=Y.*G1(:,k);
    G1(:,k)=ifft(G1(:,k));
end
%G2,G4物理空间
G2(:,1)=2*A.*conj(B);
G2(:,2)=2*A.*conj(A)+2*B.*conj(B);
G2(:,3)=2*conj(A).*B;
G4(:,1)=A.^2;
G4(:,2)=2*A.*B;
G4(:,3)=B.^2;

G7(:,1)=2*(A.*conj(F1)+conj(B).*E1);
G7(:,2)=2*(A.*conj(E1)+B.*conj(F1)+conj(B).*F1+conj(A).*E1);
G7(:,3)=2*(B.*conj(E1)+conj(A).*F1);


G10(:,1)=2*(A.*E1);
G10(:,2)=2*(A.*F1+B.*E1);
G10(:,3)=2*(B.*F1);
G12(:,1)=A;
G12(:,2)=B;

    %nonlinear term
    nontf2=zeros(M,1);
    nontg2=nontf2;
   for k=1:3
       for kk=1:4%二阶只有两项(第一行对)
 
           nontf2=nontf2+f5a(k,kk)*fft(G2(:,k).*G1(:,kk))+f6a(k,kk)*fft(G4(:,k).*conj(G1(:,kk)))...
              +fs5a(k,kk)*fft(G7(:,k).*G1(:,kk))+fs6a(k,kk)*fft(G10(:,k).*conj(G1(:,kk)))...
               +fws5a(k,kk)*i*betag.*fft(G2(:,k).*G1(:,kk))+fws6a(k,kk)*i*betag.*fft(G4(:,k).*conj(G1(:,kk)));
              
           nontg2=nontg2+f7a(k,kk)*fft(G2(:,k).*G1(:,kk))+f8a(k,kk)*fft(G4(:,k).*conj(G1(:,kk)))...
               +fs7a(k,kk)*fft(G7(:,k).*G1(:,kk))+fs8a(k,kk)*fft(G10(:,k).*conj(G1(:,kk)))...
               +fws7a(k,kk)*i*betag.*fft(G2(:,k).*G1(:,kk))+fws8a(k,kk)*i*betag.*fft(G4(:,k).*conj(G1(:,kk)));
              
       end
   end
   %以下也是三阶的
   for k=1:4
       for kk=1:4
           for kkk=1:2

           nontf2=nontf2+f10a(k,kk,kkk)*fft(2*G1(:,k).*conj(G1(:,kk)).*G12(:,kkk))...
           +f12a(k,kk,kkk)*fft(G1(:,k).*G1(:,kk).*conj(G12(:,2-kkk+1)));
     nontg2=nontg2+g10a(k,kk,kkk)*fft(2*G1(:,k).*conj(G1(:,kk)).*G12(:,kkk))...
           +g12a(k,kk,kkk)*fft(G1(:,k).*G1(:,kk).*conj(G12(:,2-kkk+1)));
       end
       end
   end

u0=u0+Y.*(nontf1+nontf2);
u1=u1+z*(nontg1+nontg2);
% third order
nontf3=zeros(M,1);
nontg3=nontf3;
for k=1:3
    for kk=1:3
        for kkk=1:4
           
            
            nontf3=nontf3+f16a(k,kk,kkk)*fft(G2(:,k).*ifft(Y.*fft(G2(:,kk).*G1(:,kkk))))+f18a(k,kk,kkk)*fft(G2(:,k).*ifft(Y.*fft(G4(:,kk).*conj(G1(:,kkk)))))...
                +f17a(k,kk,kkk)*fft(G4(:,k).*conj(ifft(Y.*fft(G2(:,kk).*G1(:,kkk)))))+f19a(k,kk,kkk)*fft(G4(:,k).*conj(ifft(Y.*fft(G4(:,kk).*conj(G1(:,kkk))))));
            nontg3=nontg3+f20a(k,kk,kkk)*fft(G2(:,k).*ifft(Y.*fft(G2(:,kk).*G1(:,kkk))))+f22a(k,kk,kkk)*fft(G2(:,k).*ifft(Y.*fft(G4(:,kk).*conj(G1(:,kkk)))))...
                +f21a(k,kk,kkk)*fft(G4(:,k).*conj(ifft(Y.*fft(G2(:,kk).*G1(:,kkk)))))+f23a(k,kk,kkk)*fft(G4(:,k).*conj(ifft(Y.*fft(G4(:,kk).*conj(G1(:,kkk))))));
        end
    end
end
for k=1:3
    for kk=1:4
       
        
        nontf3=nontf3+f24a(k,kk)*fft(G2(:,k).*ifft(Y.*G56(:,kk)))+f25a(k,kk)*fft(G4(:,k).*conj(ifft(Y.*G56(:,kk))))...
            +f28a(k,kk)*fft(G2(:,k).*ifft(i*betag.*fft(G1(:,kk))))+f29a(k,kk)*fft(G4(:,k).*conj(ifft(i*betag.*fft(G1(:,kk)))));
        nontg3=nontg3+f26a(k,kk)*fft(G2(:,k).*ifft(Y.*G56(:,kk)))+f27a(k,kk)*fft(G4(:,k).*conj(ifft(Y.*G56(:,kk))))...
            +f30a(k,kk)*fft(G2(:,k).*ifft(i*betag.*fft(G1(:,kk))))+f31a(k,kk)*fft(G4(:,k).*conj(ifft(i*betag.*fft(G1(:,kk)))));
    end
end

u0=u0+Y.*(nontf3);
u1=u1+z*(nontg3);
if mod(tk,100)==0%not save the intial value
    FU(:,(tk)/100,fe)=ifft(u0);
end
%u(:,tk)=ifft(u0);%为了像包老师一样随着时间验证
%
end

    u0=ifft(u0);
%% time end
%error analysis
% FU=u0;
% %space error
% if fN>1
%     Fe(fN-1,fe)=0;
%     for l=1:2:M
%          Fe(fN-1,fe)=Fe(fN-1,fe)+sqrt(h)*sqrt(conj(FG((l-1)/2+1)-FU(l))*(FG((l-1)/2+1)-FU(l)));
%     end
% end
% FG=u0;


%FU(:,fN,fe)=u0;
% t(fN)=tau;
% if fN>1
%     Fe(fN-1,fe)=sqrt(h)*norm(FU(:,fN-1)-FU(:,fN),2);
%      tt(fN-1)=t(fN-1)-t(fN);
% end
end%fN
% plot
disp('ended');
% plot(log2(tt),log2(Fe(:,fe)),'b-s',log2(tt),log2(tt),'k',log2(tt),log2(tt),'k-*');
% hold on

% x=linspace(a,b,M);
% plot(x,real(u0))
% y=linspace(0,T,N+1);
% [X,Y]=meshgrid(y,x);
% figure
% mesh(X,Y,real(u0))
% time periodic
% at=linspace(0,T,N+1);
% figure;
% subplot(2,3,1);
% plot(at,real(UU(1,:)))
% NFFT = 2^nextpow2(N);
% subplot(2,3,4);
% att= N/2*linspace(0,T,NFFT/2+1);
% cc=fft(UU(1,:),NFFT)/N;
% plot(att,abs(cc(1:NFFT/2+1)))
% 
% subplot(2,3,2);
% plot(at,real(UU(M/2,:)))
% subplot(2,3,5);
% cc=fft(UU(M/2,:),NFFT)/N;
% plot(att,abs(cc(1:NFFT/2+1)))
% 
% subplot(2,3,3);
% plot(at,real(UU(M/4*3,:)))
% subplot(2,3,6);
% cc=fft(UU(M/4*3,:),NFFT)/N;
% plot(att,abs(cc(1:NFFT/2+1)))


%FU(:,fe)=u0;
%u(:,fe)=u0;
end%fe
%toc
uw=FU;
end
%% 














