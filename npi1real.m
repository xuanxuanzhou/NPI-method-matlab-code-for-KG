%% handbook
%the program is for the real-value quadratic 
%Klein-Gordon equation in the non-relativistic limit regime,
%first-order Nested Picard method 
%with Fourier pseudospectral method,
%function: function [us,t]=npi1real(T)
%input: T is T_max
%output: us is the numerical solution 
%(the index of column is fe, the index of row is fN, t is the same as us),
%t is the cpu time 
%% function
clc;
clear;
syms s w epsilo real;
%first-order used funciton
p{1}=symfun(epsilo^2/(3*1i)*exp(-1i*s/epsilo^2)*(exp(3*1i*s/epsilo^2)-1),[w,epsilo]);
%symfun(exp(-1i*s/epsilo^2)*int(exp(3*1i*w/epsilo^2),w,0,s),[w,epsilo]);
p{2}=symfun(epsilo^2/(1*1i)*exp(-1i*s/epsilo^2)*(exp(1*1i*s/epsilo^2)-1),[w,epsilo]);
%symfun(exp(-1i*s/epsilo^2)*int(exp(1*1i*w/epsilo^2),w,0,s),[w,epsilo]);
p{3}=symfun(epsilo^2/(-1*1i)*exp(-1i*s/epsilo^2)*(exp(-1*1i*s/epsilo^2)-1),[w,epsilo]);
%symfun(exp(-1i*s/epsilo^2)*int(exp(-1*1i*w/epsilo^2),w,0,s),[w,epsilo]);
for k=1:3
    p{k}=subs(p{k},s,w);
end
%function [us,t]=npi1real(T)
tic
%load('symp.mat');
T=1;                                                                    %the T_max
a=-128;b=128;                                                          %the domain
lada=1;                                                               %the coefficient of nonlinear term
%ht=waitbar(0,'Please wait...');                              %the progress bar
%% preparation
for fe=1:2:5                                                          %control the value of epsilon
    epsilon=1/2^(2*(fe-1));
    for  fN=1:7                                                  %control the value of M or N
        disp(fN);
        %M=2^3*2^(fN-1);
        M=128*32;%2^6;                                                        %the spatial discretization
        h=(b-a)/M;
        %N=1000;                                                      %the temporal discretization
        N=100*2^(fN-1);
        tau=T/N; 
        ps=zeros(3,1);
        for k=1:3                                                      %integral coefficient
            ps(k)=double(p{k}(tau,epsilon));
        end
        % initialization
        u0=zeros(M,1); u1=u0;
        for k=1:M
            %the first set of initial conditions
%             u0(k,1)=1/2*(cos(3*(a+(k-1)*h))^2*sin(2*(a+(k-1)*h)))/(2-cos(a+(k-1)*h));
%             u1(k,1)=1/(2*epsilon^2)*(cos(2*(a+(k-1)*h))*sin(a+(k-1)*h))/(2-cos(a+(k-1)*h));
        u0(k,1)=exp(-(a+(k-1)*h)^2)/sqrt(pi);
        u1(k,1)=1/epsilon^2*1/2*sech((a+(k-1)*h)^2)*sin(a+(k-1)*h);
%         u0(k,1)=sech(a+(k-1)*h)^2;
%         u1(k,1)=1/epsilon^2*-sqrt(3)*sech(a+(k-1)*h)^2*tanh(a+(k-1)*h);
        end
        % sptial matrix
        mul=zeros(M,1); betal=mul;
        for l=1:M                                                     %notice the solver of 0->0(the first method for l: 0~M/2-1<=>1~M/2, -M/2~-1<=>M/2+1~M)
            if l>M/2
                mul(l)=2*pi*(l-M-1)/(b-a);
            else
                mul(l)=2*pi*(l-1)/(b-a);
            end
            betal(l,1)=sqrt(1+epsilon^2*mul(l)^2)/(epsilon^2);
        end
        u0=fft(u0); u1=fft(u1);
        A=1i*lada./(2*epsilon^2*betal); 
        up=1/2*(u0+1i./betal.*u1);
        %% computing
        for tk=1:N                                                    %notice the tk can not be used again below
            %delta^{n,1}_{+} delta_n1p
            ups=ifft(up);
            F_1=A.*fft(conj(ups).^2); F_2=2*A.*fft(ups.*conj(ups)); F_3=A.*fft(ups.^2);
            
            delta_n1p=ps(1)*F_1 + ps(2)*F_2+ ps(3)*F_3;
            
            ups=exp(-1i*tau*betal).*up+ delta_n1p;
            
            up=ups;
        end
        nonerr1(fN,fe)=norm(delta_n1p,2)/tau;
        nonerr2(fN,fe)=norm(delta_n1p,'inf')/tau;
        %waitbar(fN*fe/21,ht);                                       %show the progress bar
        us{fN,fe}=2*real(ifft(up));                                              %save the numerical solution
        if fN>1
            err(fN-1,fe)=sqrt(h)*norm(us{fN,fe}-us{fN-1,fe},2);
            ts(fN-1)=tau;
        end
    end%fN
    %us(:,fe)=up;
      t(fN,fe)=toc;                                                         %save the cpu time
end%fe
%close(ht);                                                      %close the progress bar
% plot(log2(ts),log2(err(:,fe)),'r',log2(ts),1*log2(ts))
% x=linspace(-pi,pi-2*pi/M,M);
% plot(x,up)
%end%function



