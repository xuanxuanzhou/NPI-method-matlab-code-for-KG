%% handbook
%the program is for the real-value quadratic
%Klein-Gordon equation in the non-relativistic limit regime,
%third-order Nested Picard method
%with Fourier pseudospectral method,
%function: function [us,t]=npi3real(flag)
%input: the explanation of flag can be found below
%output: us is the numerical solution
%(the index of column is fe, the index of row is fN, t is the same as us),
%t is the cpu time
%% function
function [us,t]=npi3real(flag)
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
            p2s(k)=double(p2{k}(tau,epsilon));
            q2s(k)=double(q2{k}(tau,epsilon));
            q2cs(k)=double(q2c{k}(tau,epsilon));
            q0s(k)=double(q0{k}(tau,epsilon));
            q0cs(k)=double(q0c{k}(tau,epsilon));
            
            q21s(k)=double(q21{k}(tau,epsilon));
            q2c1s(k)=double(q2c1{k}(tau,epsilon));
            q01s(k)=double(q01{k}(tau,epsilon));
            q0c1s(k)=double(q0c1{k}(tau,epsilon));
        end
        for k=1:3
            for kk=1:3
                rs(k,kk)=double(r{k,kk}(tau,epsilon));
                rs(k+3,kk)=double(r{k+3,kk}(tau,epsilon));
                rs(k,kk+3)=double(r{k,kk+3}(tau,epsilon));
                rs(k+3,kk+3)=double(r{k+3,kk+3}(tau,epsilon));
            end
        end
        
        for k=1:6
            m2s(k)=double(m2{k}(tau,epsilon));
            m2cs(k)=double(m2c{k}(tau,epsilon));
            m0s(k)=double(m0{k}(tau,epsilon));
            m0cs(k)=double(m0c{k}(tau,epsilon));
            r2s(k)=double(r2{k}(tau,epsilon));
            r2cs(k)=double(r2c{k}(tau,epsilon));
            r0s(k)=double(r0{k}(tau,epsilon));
            r0cs(k)=double(r0c{k}(tau,epsilon));
        end
        %% computing
        for tk=1:N %notice the tk can not be used again below
            delta_n1p=zeros(M,1); delta_n2p=zeros(M,1); delta_n3p=zeros(M,1);
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
            %ŒÔ¿Ìø’º‰
            for k=1:3
                F1(:,k)=ifft(F1(:,k));
            end
            for k=1:6
                E2(:,k)=ifft(E2(:,k));
            end
            %third
            B2up=ifft(B.^2.*fft(conj(up)));
            F2(:,1)=A.*fft(Bup.^2);
            F2(:,2)=2*A.*fft(abs(Bup).^2);
            F2(:,3)=A.*fft(conj(Bup).^2);
            
            F3(:,1)=A.*fft(conj(up).*B2up);
            F3(:,2)=A.*fft(conj(up).*conj(B2up) + up.*B2up);
            F3(:,3)=A.*fft(up.*conj(B2up));
            for k=1:3
                E3(:,k)=2*Bup.*F(:,k);
                E3(:,k)=fft(E3(:,k));
                H2(:,k)=2*conj(up).*F1(:,k);
                H2(:,k)=fft(H2(:,k));
                G2(:,k)=2*conj(up).*ifft(A.*fft(E2(:,k)));
                G2(:,k)=fft(G2(:,k));
            end
            for k=4:6
                E3(:,k)=2*Bup.*conj(F(:,k-3));
                E3(:,k)=fft(E3(:,k));
                H2(:,k)=2*conj(up).*conj(F1(:,k-3));
                H2(:,k)=fft(H2(:,k));
                G2(:,k)=2*conj(up).*conj(ifft(A.*fft(E2(:,k-3))));
                G2(:,k)=fft(G2(:,k));
            end
            for k=7:9
                H2(:,k)=2*conj(up).*ifft(B.*fft(F(:,k-6)));
                H2(:,k)=fft(H2(:,k));
                G2(:,k)=2*conj(up).*ifft(A.*fft(E2(:,k-3)));
                G2(:,k)=fft(G2(:,k));
                
            end
            for k=10:12
                H2(:,k)=2*conj(up).*conj(ifft(B.*fft(F(:,k-9))));
                H2(:,k)=fft(H2(:,k));
                G2(:,k)=2*conj(up).*conj(ifft(A.*fft(E2(:,k-6))));
                G2(:,k)=fft(G2(:,k));
            end
            
            for k=1:3
                for kk=1:3
                    G(:,k,kk)=F(:,k).*F(:,kk);
                    G(:,k,kk)=fft(G(:,k,kk));
%                     G(:,k+3,kk+3)=conj(F(:,k)).*F(:,kk);
%                     G(:,k+3,kk+3)=fft(G(:,k+3,kk+3));
                    G(:,k+3,kk)=conj(F(:,k)).*F(:,kk);
                    G(:,k+3,kk)=fft(G(:,k+3,kk));
                    G(:,k,kk+3)=F(:,k).*conj(F(:,kk));
                    G(:,k,kk+3)=fft(G(:,k,kk+3));
                    G(:,k+3,kk+3)=conj(F(:,k)).*conj(F(:,kk));
                    G(:,k+3,kk+3)=fft(G(:,k+3,kk+3));
                end
            end
            
            for k=1:3
                F(:,k)=fft(F(:,k));
            end
            for k=1:3
                F1(:,k)=fft(F1(:,k));
            end
            for k=1:6
                E2(:,k)=fft(E2(:,k));
            end
            for k=1:3
                delta_n3p=delta_n3p+ ps(k)*tau^2*1/2*B.^2.*F(:,k) + p2s(k)*1/2*B.^2.*F(:,k) -tau*p1s(k)*B.^2.*F(:,k)...
                    -p1s(k)*tau*B.*F1(:,k)...
                    -q2s(k)*tau*A.*B.*E2(:,k)- q2cs(k)*tau*A.*B.*E2(:,k+3)- q0s(k)*tau*A.*B.*fft(conj(ifft(E2(:,k+3))))- q0cs(k)*tau*A.*B.*fft(conj(ifft(E2(:,k))))...
                    +q21s(k)*A.*B.*E2(:,k)+ q2c1s(k)*A.*B.*E2(:,k+3)+ q01s(k)*A.*B.*fft(conj(ifft(E2(:,k+3))))+ q0c1s(k)*A.*B.*fft(conj(ifft(E2(:,k))))...
                    +p2s(k)*( B.*F1(:,k) +F2(:,k) +F3(:,k)  )...
                    +q21s(k)*A.*E3(:,k)+ q2c1s(k)*A.*E3(:,k+3)+ q01s(k)*A.*fft(conj(ifft(E3(:,k+3))))+ q0c1s(k)*A.*fft(conj(ifft(E3(:,k))))...
                    +rs(k,kk)*A.*G(:,k,kk) + rs(k+3,kk)*A.*G(:,k+3,kk)+ rs(k,kk+3)*A.*G(:,k,kk+3)+ rs(k+3,kk+3)*A.*G(:,k+3,kk+3)...
                    +m2s(k)*A.*H2(:,k)+ m2cs(k)*A.*H2(:,k+3)+ m0s(k)*A.*fft(conj(ifft(H2(:,k+3))))+ m0cs(k)*A.*fft(conj(ifft(H2(:,k))))...
                    +m2s(k+3)*A.*H2(:,k+6)+ m2cs(k+3)*A.*H2(:,k+9)+ m0s(k+3)*A.*fft(conj(ifft(H2(:,k+9))))+ m0cs(k+3)*A.*fft(conj(ifft(H2(:,k+6))))...
                    +r2s(k)*A.*G2(:,k)+ r2cs(k)*A.*G2(:,k+3)+ r0s(k)*A.*fft(conj(ifft(G2(:,k+3))))+ r0cs(k)*A.*fft(conj(ifft(G2(:,k))))...
                    +r2s(k+3)*A.*G2(:,k+6)+ r2cs(k+3)*A.*G2(:,k+9)+ r0s(k+3)*A.*fft(conj(ifft(G2(:,k+9))))+ r0cs(k+3)*A.*fft(conj(ifft(G2(:,k+6))));
            end
            
            up=fft(up);
            up=exp(-1i*tau*betal).*up+ delta_n1p+ delta_n2p+ delta_n3p;
        end
        
        up=2*real(ifft(up));
        %% end
        us{fN,fe}=up;
        %waitbar(fe*fN/14,ht);
    end%fN
    t(fN,fe)=toc;
end%fe
%close(ht);

end%function



