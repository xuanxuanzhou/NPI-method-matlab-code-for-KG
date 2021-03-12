%% handbook
%the program is for the real-value quadratic 
%Klein-Gordon equation in the non-relativistic limit regime,
%second-order Nested Picard method 
%with Fourier pseudospectral method,
%function: function [us,t]=npi2real(T)
%input: T is T_max
%output: us is the numerical solution 
%(the index of column is fe, the index of row is fN, t is the same as us),
%t is the cpu time 
%% function
%function [err,us,t]=npi2real(T)
tic
% load('symp.mat');
T=1;                                                                    %the T_max
a=-32;b=32;                                                          %the domain
lada=1;                                                               %the coefficient of nonlinear term
ht=waitbar(0,'Please wait...');                              %the progress bar
%% preparation
for fe=1:4                                                       %control the value of epsilon
    epsilon=1/10^(1*(fe));
    for  fN=1:1                                                       %control the value of M or N
        %M=2^3*2^(fN-1);
        M=32*32;                                                        %the spatial discretization
        h=(b-a)/M;
        N=100000;                                                    %the temporal discretization
        %N=100*2^(fN-1);
        tau=T/N; tt=tau;                                           %tt is used by sin(tt*D)/tt
        % 符号积分系数
        ps=zeros(3,1);p1s=ps;q2s=ps;q2ls=ps;q0s=ps;q0ls=ps;
        for k=1:3
            ps(k)=double(p{k}(tau,epsilon));
            p1s(k)=double(p1{k}(tau,epsilon));
            q2s(k)=double(q2{k}(tau,epsilon));
            q2ls(k)=double(q2l{k}(tau,epsilon));
            q0s(k)=double(q0{k}(tau,epsilon));
            q0ls(k)=double(q0l{k}(tau,epsilon));
        end
        
        %initialization
        u0=zeros(M,1); u1=u0;
        for k=1:M
            %the first set of initial conditions
%             u0(k,1)=1/2*(cos(3*(a+(k-1)*h))^2*sin(2*(a+(k-1)*h)))/(2-cos(a+(k-1)*h));
%             u1(k,1)=1/(2*epsilon^2)*(cos(2*(a+(k-1)*h))*sin(a+(k-1)*h))/(2-cos(a+(k-1)*h));
%         u0(k,1)=sech(a+(k-1)*h)^2;
%         u1(k,1)=1/epsilon^2*-sqrt(3)*sech(a+(k-1)*h)^2*tanh(a+(k-1)*h);
         u0(k,1)=exp(-(a+(k-1)*h)^2)/sqrt(pi);
        u1(k,1)=1/epsilon^2*1/2*sech((a+(k-1)*h)^2)*sin(a+(k-1)*h);
        end
        % sptial matrix
        mul=zeros(M,1); betal=mul;
        for l=1:M %notice the solver of 0->0(the first method for l: 0~M/2-1<=>1~M/2, -M/2~-1<=>M/2+1~M)
            if l>M/2
                mul(l)=2*pi*(l-M-1)/(b-a);
            else
                mul(l)=2*pi*(l-1)/(b-a);
            end
            betal(l,1)=sqrt(1+epsilon^2*mul(l)^2)/(epsilon^2);
        end
        %computing preparation
        u0=fft(u0); u1=fft(u1);
        A=1i*lada./(2*epsilon^2*betal);
        %B=i*(betal-1/epsilon^2);
        B=1i*sin(tt*(betal-1/epsilon^2))/tt;
        up=1/2*(u0+1i./betal.*u1);
        %% computing
        for tk=1:N %notice the tk can not be used again below
            %delta^{n,1}_{+} delta_n1p
            up=ifft(up);
            F_1=A.*fft(conj(up).^2);
            F_2=2*A.*fft(abs(up).^2);
            F_3=A.*fft(up.^2);
            
            Bup=ifft(B.*fft(conj(up)));
            F_11=2*A.*fft(conj(up).*Bup);
            F_21=2*A.*fft(conj(up).*conj(Bup) + up.*Bup);
            F_31=2*A.*fft(up.*conj(Bup));
            
            E_21=2*fft(conj(up).*ifft(F_1));
            E_22=2*fft(conj(up).*ifft(F_2));
            E_23=2*fft(conj(up).*ifft(F_3));
            E_24=2*fft(conj(up).*conj(ifft(F_1)));
            E_25=2*fft(conj(up).*conj(ifft(F_2)));
            E_26=2*fft(conj(up).*conj(ifft(F_3)));
            
            delta_n1p=ps(1)*F_1 + ps(2)*F_2+ ps(3)*F_3;
            delta_n2p=p1s(1)*F_11 + p1s(2)*F_21+ p1s(3)*F_31 + p1s(1)*B.*F_1 + p1s(2)*B.*F_2+ p1s(3)*B.*F_3...
                -tau*ps(1)*B.*F_1 - tau*ps(2)*B.*F_2- tau*ps(3)*B.*F_3...
                +q2s(1)*A.*E_21+ q2s(2)*A.*E_22 + q2s(3)*A.*E_23 +q2ls(1)*A.*E_24+ q2ls(2)*A.*E_25 + q2ls(3)*A.*E_26...
                 +q0s(1)*A.*fft(conj(ifft(E_24)))+ q0s(2)*A.*fft(conj(ifft(E_25))) + q0s(3)*A.*fft(conj(ifft(E_26)))...
                 +q0ls(1)*A.*fft(conj(ifft(E_21)))+ q0ls(2)*A.*fft(conj(ifft(E_22))) + q0ls(3)*A.*fft(conj(ifft(E_23)));
            up=fft(up);
            % first order
            for l=1:M
                up(l)=exp(-1i*tau*betal(l))*up(l)+ delta_n1p(l)+ delta_n2p(l);
            end
             if mod(tk,1000)==0
        (tk)/1000
        FU(:,(tk)/1000,fe)=2*real(ifft(up));
           end
            
        end
        
        %up=2*real(ifft(up));
        %% end
        %us2{fN,fe}=up;%save the numerical solution
%         if fN>1
%             fN
%             err(fN-1,fe)=sqrt(h)*norm(us2{fN,fe}-us2{fN-1,fe},'fro')
%         end
     %show the progress bar
     %waitbar(fe*fN/36,ht);
    end%fN
    %u(:,fe)=up;

end%fe
%close the progress bar
%close(ht);
%t(fN,fe)=toc;
% x=linspace(-pi,pi-2*pi/M,M);
% plot(x,up)
%end%function



