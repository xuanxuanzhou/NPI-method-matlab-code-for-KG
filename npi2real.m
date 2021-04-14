%% handbook
%the program is for the real-value quadratic
%Klein-Gordon equation in the non-relativistic limit regime,
%second-order Nested Picard method
%with Fourier pseudospectral method,
%function: function [us,t]=npi2real(flag)
%input: the explanation of flag can be found below
%output: us is the numerical solution
%(the index of column is fe, the index of row is fN, t is the same as us),
%t is the cpu time
%% function
function [us,t]=npi2real(flag)
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
        M=32*32;
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
                % flag1=2,  [a,b]=[-32,32]   M is changed with a,b
            case 2
                a=-32;b=32;
                h=(b-a)/M;
                x=linspace(a,b-h,M);
                u0(:,1)=3*sin(x)./(exp(x.^2/2)+exp(-x.^2/2) );
                u1(:,1)=1/epsilon^2*2*exp(-x.^2)/sqrt(pi);
                % flag1=3,  [a,b]=[-128,128]   M is changed with a,b
            case 3
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
        A=1i*lada./(2*epsilon^2*betal);
        B=1i*sin(tt*(betal-1/epsilon^2))/tt;
        up=1/2*(u0+1i./betal.*u1);
        %upload the symbol integration
        load('symbol_intereal');
        for k=1:3
            ps(k)=double(p{k}(tau,epsilon));
            p1s(k)=double(p1{k}(tau,epsilon));
            q2s(k)=double(q2{k}(tau,epsilon));
            q2cs(k)=double(q2c{k}(tau,epsilon));
            q0s(k)=double(q0{k}(tau,epsilon));
            q0cs(k)=double(q0c{k}(tau,epsilon));
        end
        %% computing
        for tk=1:N %notice the tk can not be used again below
            %delta^{n,1}_{+} delta_n1p
            %             up=ifft(up);
            %             F_1=A.*fft(conj(up).^2);
            %             F_2=2*A.*fft(abs(up).^2);
            %             F_3=A.*fft(up.^2);
            %
            %             Bup=ifft(B.*fft(conj(up)));
            %             F_11=2*A.*fft(conj(up).*Bup);
            %             F_21=2*A.*fft(conj(up).*conj(Bup) + up.*Bup);
            %             F_31=2*A.*fft(up.*conj(Bup));
            %
            %             E_21=2*fft(conj(up).*ifft(F_1));
            %             E_22=2*fft(conj(up).*ifft(F_2));
            %             E_23=2*fft(conj(up).*ifft(F_3));
            %             E_24=2*fft(conj(up).*conj(ifft(F_1)));
            %             E_25=2*fft(conj(up).*conj(ifft(F_2)));
            %             E_26=2*fft(conj(up).*conj(ifft(F_3)));
            %
            %             delta_n1p=ps(1)*F_1 + ps(2)*F_2+ ps(3)*F_3;
            %             delta_n2p=p1s(1)*F_11 + p1s(2)*F_21+ p1s(3)*F_31 + p1s(1)*B.*F_1 + p1s(2)*B.*F_2+ p1s(3)*B.*F_3...
            %                 -tau*ps(1)*B.*F_1 - tau*ps(2)*B.*F_2- tau*ps(3)*B.*F_3...
            %                 +q2s(1)*A.*E_21+ q2s(2)*A.*E_22 + q2s(3)*A.*E_23 +q2ls(1)*A.*E_24+ q2ls(2)*A.*E_25 + q2ls(3)*A.*E_26...
            %                  +q0s(1)*A.*fft(conj(ifft(E_24)))+ q0s(2)*A.*fft(conj(ifft(E_25))) + q0s(3)*A.*fft(conj(ifft(E_26)))...
            %                  +q0ls(1)*A.*fft(conj(ifft(E_21)))+ q0ls(2)*A.*fft(conj(ifft(E_22))) + q0ls(3)*A.*fft(conj(ifft(E_23)));
            %             up=fft(up);
            %
            %             for l=1:M
            %                 up(l)=exp(-1i*tau*betal(l))*up(l)+ delta_n1p(l)+ delta_n2p(l);
            %             end
            delta_n1p=zeros(M,1); delta_n2p=zeros(M,1);
            %delta^{n,1}_{+} delta_n1p
            up=ifft(up);
            F(:,1)=A.*fft(conj(up).^2);
            F(:,2)=2*A.*fft(abs(up).^2);
            F(:,3)=A.*fft(up.^2);
            delta_n1p=ps(1)*F(:,1) + ps(2)*F(:,2)+ ps(3)*F(:,3);
            for k=1:3
                F(:,k)=ifft(F(:,k));
            end
            %second
            Bup=ifft(B.*fft(conj(up)));
            
            F1(:,1)=2*A.*fft(conj(up).*Bup);
            F1(:,2)=2*A.*fft(conj(up).*conj(Bup) + up.*Bup);
            F1(:,3)=2*A.*fft(up.*conj(Bup));
            E2(:,1)=2*fft(conj(up).*F(:,1));
            E2(:,2)=2*fft(conj(up).*F(:,2));
            E2(:,3)=2*fft(conj(up).*F(:,3));
            E2(:,4)=2*fft(conj(up).*conj(F(:,1)));
            E2(:,5)=2*fft(conj(up).*conj(F(:,2)));
            E2(:,6)=2*fft(conj(up).*conj(F(:,3)));
            for k=1:3
                delta_n2p=delta_n2p+p1s(k)*F1(:,k) +  p1s(k)*B.*fft(F(:,k)) -tau*ps(k)*B.*fft(F(:,k)) ...
                    +q2s(k)*A.*E2(:,k)+ q2cs(k)*A.*E2(:,k+3)+q0s(k)*A.*fft(conj(ifft(E2(:,k+3))))+q0cs(k)*A.*fft(conj(ifft(E2(:,k))));
            end
            up=fft(up);
            up=exp(-1i*tau*betal).*up+ delta_n1p+ delta_n2p;
        end
        us{fN,fe}=2*real(ifft(up));%save the numerical solution
        %show the progress bar
        %waitbar(fe*fN/36,ht);
    end%fN
    t(fN,fe)=toc;                                                         %save the cpu time
end%fe
%close the progress bar
%close(ht);
end%function



