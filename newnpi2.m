%% handbook
%the program is for the complex-value |u|^2u
%Klein-Gordon equation in the non-relativistic limit regime,
%second-order Nested Picard method
%with Fourier pseudospectral method,
%function: function [us,t]=newnpi2(flag)
%input: the explanation of flag can be found below
%output: us is the numerical solution
%(the index of column is fe, the index of row is fN, t is the same as us),
%t is the cpu time
%% function
function [us,t]=newnpi2(flag)
tic
%% preperation
T=1;
%the T_max
lada=1;
%the coefficient of nonlinear term
%ht=waitbar(0,'Please wait...');
%the progress bar
%% fe
for fe=1:2
    %control the value of epsilon
    epsilon=1/2^(2*(fe-1));
    for  fN=1:7
        %control the value of M or N
        %M=2^3*2^(fN-1);
        M=32*16;
        %the spatial discretization
        
        %the temporal discretization
        N=10*2^(fN-1);
        %N=1000;
        tau=T/N; tt=tau;
        %% initial value
        flag1=flag;
        switch flag1
            %the value of flag1 can be 1,2,3,4
            % flag1=1,  [a,b]=[-pi,pi]   M is changed with a,b
            case 1
                a=-pi;b=pi;
                h=(b-a)/M;
                x=linspace(a,b-h,M);
                u0(:,1)=1/2*(cos(3*x).^2.*sin(2*x))./(2-cos(x));
                u1(:,1)=1/(2*epsilon^2)*(cos(2*x).*sin(x))./(2-cos(x));
                % flag1=2,  [a,b]=[-pi,pi]   M is changed with a,b
            case 2
                a=-pi;b=pi;
                h=(b-a)/M;
                x=linspace(a,b-h,M);
                u0(:,1)=(2+i)/sqrt(5)*cos(x);
                u1(:,1)=1/epsilon^2*((1+i)/sqrt(2)*sin(x)+1/2*cos(x));
                % flag1=3,  [a,b]=[-32,32]   M is changed with a,b
            case 3
                a=-32;b=32;
                h=(b-a)/M;
                x=linspace(a,b-h,M);
                u0(:,1)=3*sin(x)./(exp(x.^2/2)+exp(-x.^2/2) );
                u1(:,1)=1/epsilon^2*2*exp(-x.^2)/sqrt(pi);
                % flag1=4,  [a,b]=[-128,128]   M is changed with a,b
            case 4
                a=-32;b=32;
                h=(b-a)/M;
                x=linspace(a,b-h,M);
                u0(:,1)=exp(-x.^2)/sqrt(pi);
                u1(:,1)=1/epsilon^2*1/2*sech(x.^2).*sin(x);
            otherwise
                disp('wrong number!');
        end
        
        mul=zeros(M,1); betal=mul;
        for l=1:M
            %notice the solver of 0->0(the first method for l: 0~M/2-1<=>1~M/2, -M/2~-1<=>M/2+1~M)
            if l>M/2
                mul(l)=2*pi*(l-M-1)/(b-a);
            else
                mul(l)=2*pi*(l-1)/(b-a);
            end
            betal(l,1)=sqrt(1+epsilon^2*mul(l)^2)/(epsilon^2);
        end
        u0=fft(u0); u1=fft(u1);
        betag=betal-1/epsilon^2;
        B=i*sin(tt*(betag))/tt;
        A=1i*lada./(2*epsilon^2*betal);
        up=1/2*(u0+1i./betal.*u1);
        um=1/2*(u0-1i./betal.*u1);
        %upload the symbol integration
        load('symbol_inte');
        for k=1:4                                                      %integral coefficient
            ps(k)=double(p{k}(tau,epsilon));
            p1s(k)=double(p1{k}(tau,epsilon));
        end
        for j=1:6
            for k=1:2
                qs(j,k)=double(q{j,k}(tau,epsilon));
            end
        end
        %% computing
        for tk=1:N                                                    %notice the tk can not be used again below
            %delta^{n,1}_{+} delta_n1p
            ups=ifft(up);
            ums=ifft(um);
            ups1=ifft(B.*up);
            ums1=ifft(B.*um);
            %% first-order
            F_1p=ifft(A.*fft(ups.^2.*conj(ums)));
            F_1m=ifft(A.*fft(ums.^2.*conj(ups)));
            F_2p=ifft(A.*fft((2*abs(ums).^2+abs(ups).^2).*ups));
            F_2m=ifft(A.*fft((2*abs(ups).^2+abs(ums).^2).*ums));
            %% second-order
            G_1p=A.*fft(ifft(B.*fft(ups.^2.*conj(ums)))   );
            G_1m=A.*fft(ifft(B.*fft(ums.^2.*conj(ups))) );
            G_2p=A.*fft( ifft(B.*fft( (2*abs(ums).^2+ abs(ups).^2).*ups)) );
            G_2m=A.*fft( ifft(B.*fft( (2*abs(ups).^2+ abs(ums).^2).*ums)));
            
            G_3p=A.*fft(  ups.^2.*conj(ums1)- 2*ups1.*ups.*conj(ums)     );
            G_3m=-A.*fft( ums.^2.*conj(ups1)- 2*ums1.*ums.*conj(ups)     );
            G_4p=A.*fft( 2*ups.*conj(ums).*ums1- 2*( abs(ups).^2+ abs(ums).^2).*ups1-...
                ups.^2.*conj(ups1)+ 2*ups.*ums.*conj(ums1)   );
            G_4m=-A.*fft( 2*ums.*conj(ups).*ups1- 2*( abs(ums).^2+ abs(ups).^2).*ums1-...
                ums.^2.*conj(ums1)+ 2*ums.*ups.*conj(ups1)  );
            
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
            
            delta_n2p=p1s(1)*(G_2p+G_4p) + p1s(2)*(G_2m+G_4m)+ p1s(3)*(G_1p+ G_3p) + p1s(4)*(G_1m+G_3m)...
                -tau*ps(1)*B.*fft(F_2p) - tau*ps(2)*B.*fft(F_2m)- tau*ps(3)*B.*fft(F_1p) - tau*ps(4)*B.*fft(F_1m)+...
                qs(1,1)*F_21p+ qs(1,2)*F_21m +qs(2,1)*F_22p+ qs(2,2)*F_22m +qs(3,1)*F_23p+ qs(3,2)*F_23m+...
                qs(4,1)*F_24p+ qs(4,2)*F_24m +qs(5,1)*F_25p+ qs(5,2)*F_25m +qs(6,1)*F_26p+ qs(6,2)*F_26m;
            
            delta_n2m=conj(p1s(1))*(G_2m-G_4m) + conj(p1s(2))*(G_2p- G_4p)+ conj(p1s(3))*(G_1m- G_3m) + conj(p1s(4))*(G_1p-G_3p)...
                -tau*conj(ps(1))*B.*fft(F_2m) - tau*conj(ps(2))*B.*fft(F_2p)- tau*conj(ps(3))*B.*fft(F_1m) - tau*conj(ps(4))*B.*fft(F_1p)+...
                conj(qs(1,1))*F_21m+ conj(qs(1,2))*F_21p + conj(qs(2,1))*F_22m+ conj(qs(2,2))*F_22p + conj(qs(3,1))*F_23m+ conj(qs(3,2))*F_23p+...
                conj(qs(4,1))*F_24m+ conj(qs(4,2))*F_24p + conj(qs(5,1))*F_25m+ conj(qs(5,2))*F_25p + conj(qs(6,1))*F_26m+ conj(qs(6,2))*F_26p;
            
            ups=exp(-1i*tau*betal).*up+ fft(delta_n1p)+ delta_n2p;
            ums=exp(1i*tau*betal).*um+ fft(delta_n1m)+ delta_n2m;
            up=ups;
            um=ums;
        end
        % waitbar(fN*fe/21,ht);                                       %show the progress bar
        us{fN,fe}=ifft(up+um);                                              %save the numerical solution
    end%fN
    t(fN,fe)=toc;                                                         %save the cpu time
end%fe
%close(ht);      
%close the progress bar
end%function



