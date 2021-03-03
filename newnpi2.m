%% handbook
%the program is for the complex-value |u|^2u 
%Klein-Gordon equation in the non-relativistic limit regime,
%first-order Nested Picard method 
%with Fourier pseudospectral method,
%function: function [us,t]=newnpi1(T)
%input: T is T_max
%output: us is the numerical solution 
%(the index of column is fe, the index of row is fN, t is the same as us),
%t is the cpu time 
%% function

%function [us,t]=newnpi2(T)
tic
%load('symp.mat');
T=1;                                                                    %the T_max
a=-128;b=128;                                                          %the domain
lada=1;                                                               %the coefficient of nonlinear term
%ht=waitbar(0,'Please wait...');                              %the progress bar
%% preparation
for fe=6:2:6                                                          %control the value of epsilon
    epsilon=1/2^(2*(fe-1));
    for  fN=1:9                                                    %control the value of M or N
        disp(fN);
        %M=2^3*2^(fN-1);
        M=128*32;%2^6;                                                        %the spatial discretization
        h=(b-a)/M;
        %N=1000;                                                      %the temporal discretization
        N=10*2^(fN-1);
        tau=T/N; tt=tau;
        ps=zeros(4,1);
        for k=1:4                                                      %integral coefficient
            ps(k)=double(p{k}(tau,epsilon));
            p1s(k)=double(p1{k}(tau,epsilon));
        end
        for k=1:2
            q1s(k)=double(q1{k}(tau,epsilon));
             q2s(k)=double(q2{k}(tau,epsilon));
              q3s(k)=double(q3{k}(tau,epsilon));
               q4s(k)=double(q4{k}(tau,epsilon));
                q5s(k)=double(q5{k}(tau,epsilon));
                 q6s(k)=double(q6{k}(tau,epsilon));
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
        betag=betal-1/epsilon^2;
        u0=fft(u0); u1=fft(u1);
        A=1i*lada./(2*epsilon^2*betal); 
        B=i*sin(tt*(betag))/tt;
        up=1/2*(u0+1i./betal.*u1);
        um=1/2*(u0-1i./betal.*u1);
        %% computing
        for tk=1:N                                                    %notice the tk can not be used again below
            %delta^{n,1}_{+} delta_n1p
            ups=ifft(up);
            ums=ifft(um);
            ups1=ifft(B.*up);
            ums1=ifft(B.*um);
            
            F_1p=ifft(A.*fft(ups.^2.*conj(ums))); 
            F_1m=ifft(A.*fft(ums.^2.*conj(ups)));
            F_2p=ifft(A.*fft((2*abs(ums).^2+abs(ups).^2).*ups));
            F_2m=ifft(A.*fft((2*abs(ups).^2+abs(ums).^2).*ums));
            
            G_1p=A.*fft(ifft(B.*fft(ups.^2.*conj(ums)))+ ( ups.^2.*conj(ums1)- 2*ups1.*ups.*conj(ums)   )   );
            G_1m=A.*fft(ifft(B.*fft(ums.^2.*conj(ups)))- ( ums.^2.*conj(ups1)- 2*ums1.*ums.*conj(ups)   )   );
            G_2p=A.*fft( ifft(B.*fft( (2*abs(ums).^2+ abs(ups).^2).*ups))+ ( 2*ups.*conj(ums).*ums1- 2*( abs(ups).^2+ abs(ums).^2).*ups1-... 
                      ups.^2.*conj(ups1)+ 2*ups.*ums.*conj(ums1) )  );
            G_2m=A.*fft( ifft(B.*fft( (2*abs(ups).^2+ abs(ums).^2).*ums))- ( 2*ums.*conj(ups).*ups1- 2*( abs(ums).^2+ abs(ups).^2).*ums1-... 
                      ums.^2.*conj(ums1)+ 2*ums.*ups.*conj(ups1) )  );
                  
             F_21p=A.*fft( 2*(abs(ups).^2+abs(ums).^2).*F_1m- 2*ups.*ums.*conj(F_1p));
             F_21m=A.*fft( 2*(abs(ups).^2+abs(ums).^2).*F_1p- 2*ups.*ums.*conj(F_1m));
             F_22p=A.*fft( 2*(abs(ups).^2+abs(ums).^2).*F_2m- 2*ups.*ums.*conj(F_2p));
             F_22m=A.*fft( 2*(abs(ups).^2+abs(ums).^2).*F_2p- 2*ups.*ums.*conj(F_2m));
             
             F_23p=A.*fft( 2*ums.*conj(ups).*F_1m- ums.^2.*conj(F_1p));
             F_23m=A.*fft( 2*ups.*conj(ums).*F_1p- ups.^2.*conj(F_1m));
             F_24p=A.*fft( 2*ums.*conj(ups).*F_2m- ums.^2.*conj(F_2p));
             F_24m=A.*fft( 2*ups.*conj(ums).*F_2p- ups.^2.*conj(F_2m));
             
             F_25p=A.*fft( 2*ups.*conj(ums).*F_1m- ups.^2.*conj(F_1p));
             F_25m=A.*fft( 2*ums.*conj(ups).*F_1p- ums.^2.*conj(F_1m));
             F_26p=A.*fft( 2*ups.*conj(ums).*F_2m- ups.^2.*conj(F_2p));
             F_26m=A.*fft( 2*ums.*conj(ups).*F_2p- ums.^2.*conj(F_2m));
             
            delta_n1p=ps(1)*F_2p + ps(2)*F_2m+ ps(3)*F_1p + ps(4)*F_1m;
            delta_n1m=-conj(ps(1))*F_2m + -conj(ps(2))*F_2p+ -conj(ps(3))*F_1m + -conj(ps(4))*F_1p;
            
            delta_n2p=p1s(1)*G_2p + p1s(2)*G_2m+ p1s(3)*G_1p + p1s(4)*G_1m...
                             -tau*ps(1)*B.*fft(F_2p) - tau*ps(2)*B.*fft(F_2m)- tau*ps(3)*B.*fft(F_1p) - tau*ps(4)*B.*fft(F_1m)+...         
                             q1s(1)*F_21p+ q1s(2)*F_21m +q2s(1)*F_22p+ q2s(2)*F_22m +q3s(1)*F_23p+ q3s(2)*F_23m+...
                             q4s(1)*F_24p+ q4s(2)*F_24m +q5s(1)*F_25p+ q5s(2)*F_25m +q6s(1)*F_26p+ q6s(2)*F_26m;
                         delta_n2m=fft(conj(ifft(delta_n2p)));
%             delta_n2m=conj(p1s(1))*G_2m + conj(p1s(2))*G_2p+ conj(p1s(3))*G_1m + conj(p1s(4))*G_1p...
%                              -tau*conj(ps(1))*B.*fft(F_2m) - tau*conj(ps(2))*B.*fft(F_2p)- tau*conj(ps(3))*B.*fft(F_1m) - tau*conj(ps(4))*B.*fft(F_1p)+...         
%                              conj(q1s(1))*F_21m+ conj(q1s(2))*F_21p + conj(q2s(1))*F_22m+ conj(q2s(2))*F_22p + conj(q3s(1))*F_23m+ conj(q3s(2))*F_23p+...
%                              conj(q4s(1))*F_24m+ conj(q4s(2))*F_24p + conj(q5s(1))*F_25m+ conj(q5s(2))*F_25p + conj(q6s(1))*F_26m+ conj(q6s(2))*F_26p;
            

            ups=exp(-1i*tau*betal).*up+ fft(delta_n1p)+ delta_n2p;
            ums=exp(1i*tau*betal).*um+ fft(delta_n1m)+ delta_n2m;
            up=ups;
            um=ums;
        end
        nonerr1(fN,fe)=norm(delta_n1p+delta_n1m,inf)/tau;
        nonerr2(fN,fe)=norm(delta_n2p+delta_n2m,'inf')/tau^2;
       % waitbar(fN*fe/21,ht);                                       %show the progress bar
        us{fN,fe}=ifft(up+um);                                              %save the numerical solution
        if fN>1
            err(fN-1,fe)=sqrt(h)*norm(us{fN,fe}-us{fN-1,fe},'fro');
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



