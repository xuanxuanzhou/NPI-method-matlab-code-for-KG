%% program handbook
%the program for the Klein-Gordon equation in the non-relativistic limit regime
%second-order Nested Picard method for real quadratic system
%共轭运算和fft,ifft不可交换
%% function
% clc;clear;
%function [us,t]=npi3real(T)
tic
T=1; a=-pi; b=pi; lada=1;%the coefficient of nonlinear term
%ht= waitbar(0,'Please wait...');
%% starting
for fe=1:6
    %epsilon=10^(-(fe-1));
    epsilon=1/2^(2*(fe-1));
    %%
    for  fN=1:1
        disp(fN);
        %M=2^3*2^(fN-1);
        M=2^8;
        %M=128*2*16;
        h=(b-a)/M;
        N=10000;
        %N=10*2^(fN-1);
        tau=T/N; tt=tau;
        % 符号积分系数
        for k=1:3
            ps(k)=double(p{k}(tau,epsilon));
            p1s(k)=double(p1{k}(tau,epsilon));
            p2s(k)=double(p2{k}(tau,epsilon));
            q2s(k)=double(q2{k}(tau,epsilon));
            q2ls(k)=double(q2l{k}(tau,epsilon));
            q0s(k)=double(q0{k}(tau,epsilon));
            q0ls(k)=double(q0l{k}(tau,epsilon));
            
            q21s(k)=double(q21{k}(tau,epsilon));
            q2l1s(k)=double(q2l1{k}(tau,epsilon));
            q01s(k)=double(q01{k}(tau,epsilon));
            q0l1s(k)=double(q0l1{k}(tau,epsilon));
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
            m2ls(k)=double(m2l{k}(tau,epsilon));
            m0s(k)=double(m0{k}(tau,epsilon));
            m0ls(k)=double(m0l{k}(tau,epsilon));
            r2s(k)=double(r2{k}(tau,epsilon));
            r2ls(k)=double(r2l{k}(tau,epsilon));
            r0s(k)=double(r0{k}(tau,epsilon));
            r0ls(k)=double(r0l{k}(tau,epsilon));
        end
        
        %initialization
        u0=zeros(M,1); u1=u0;
        for k=1:M
            %the first set of initial conditions相似文章方程一致
            u0(k,1)=1/2*(cos(3*(a+(k-1)*h))^2*sin(2*(a+(k-1)*h)))/(2-cos(a+(k-1)*h));
            u1(k,1)=1/(2*epsilon^2)*(cos(2*(a+(k-1)*h))*sin(a+(k-1)*h))/(2-cos(a+(k-1)*h));%一阶导数初值
            %exact solution
            %  u0(k,1)=3*sin(a+(k-1)*h)/(exp((a+(k-1)*h)^2/2)+exp(-(a+(k-1)*h)^2/2) );
            %         u1(k,1)=1/epsilon^2*2*exp(-(a+(k-1)*h)^2)/sqrt(pi);
            %  u0(k,1)=exp(-(a+(k-1)*h)^2)/sqrt(pi);
            %         u1(k,1)=1/epsilon^2*1/2*sech((a+(k-1)*h)^2)*sin(a+(k-1)*h);
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
        A=i*lada./(2*epsilon^2*betal);
        B=i*sin(tt*(betal-1/epsilon^2))/tt;
        up=1/2*(u0+i./betal.*u1);
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
                    +q2s(k)*A.*E2(:,k)+ q2ls(k)*A.*E2(:,k+3)+q0s(k)*A.*fft(conj(ifft(E2(:,k+3))))+q0ls(k)*A.*fft(conj(ifft(E2(:,k))));
            end
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
                    G(:,k+3,kk)=conj(F(:,k)).*F(:,kk);
                    G(:,k+3,kk)=fft(G(:,k,kk));
                    G(:,k,kk+3)=F(:,k).*conj(F(:,kk));
                    G(:,k,kk+3)=fft(G(:,k,kk));
                    G(:,k+3,kk+3)=conj(F(:,k)).*conj(F(:,kk));
                    G(:,k+3,kk+3)=fft(G(:,k,kk));
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
                delta_n3p=delta_n3p+ps(k)*tau^2*1/2*B.^2.*F(:,k) + p2s(k)*1/2*B.^2.*F(:,k) -tau*p1s(k)*B.^2.*F(:,k)...
                    -p1s(k)*tau*B.*F1(:,k)...
                    -q2s(k)*tau*A.*B.*E2(:,k)- q2ls(k)*tau*A.*B.*E2(:,k+3)- q0s(k)*tau*A.*B.*fft(conj(ifft(E2(:,k+3))))- q0ls(k)*tau*A.*B.*fft(conj(ifft(E2(:,k))))...
                    +q21s(k)*tau*A.*B.*E2(:,k)+ q2l1s(k)*tau*A.*B.*E2(:,k+3)+ q01s(k)*tau*A.*B.*fft(conj(ifft(E2(:,k+3))))+ q0l1s(k)*tau*A.*B.*fft(conj(ifft(E2(:,k))))...
                    +p2s(k)*( B.*F1(:,k) +F2(:,k) +F3(:,k)  )...
                    +q21s(k)*A.*E3(:,k)+ q2l1s(k).*A.*E3(:,k+3)+ q01s(k)*A.*fft(conj(ifft(E3(:,k+3))))+ q0l1s(k).*A.*fft(conj(ifft(E3(:,k))))...
                    +m2s(k)*A.*H2(:,k)+ m2ls(k).*A.*H2(:,k+3)+ m0s(k)*A.*fft(conj(ifft(H2(:,k+3))))+ m0ls(k).*A.*fft(conj(ifft(H2(:,k))))...
                    +m2s(k+3)*A.*H2(:,k+6)+ m2ls(k+3).*A.*H2(:,k+9)+ m0s(k+3)*A.*fft(conj(ifft(H2(:,k+9))))+ m0ls(k+3).*A.*fft(conj(ifft(H2(:,k+6))))...
                    +r2s(k)*A.*G2(:,k)+ r2ls(k).*A.*G2(:,k+3)+ r0s(k)*A.*fft(conj(ifft(G2(:,k+3))))+ r0ls(k).*A.*fft(conj(ifft(G2(:,k))))...
                    +r2s(k+3)*A.*G2(:,k+6)+ r2ls(k+3).*A.*G2(:,k+9)+ r0s(k+3)*A.*fft(conj(ifft(G2(:,k+9))))+ r0ls(k+3).*A.*fft(conj(ifft(G2(:,k+6))))...
                    +rs(k,kk)*A.*G(:,k,kk) + rs(k+3,kk)*A.*G(:,k+3,kk)+ rs(k,kk+3)*A.*G(:,k,kk+3)+ rs(k+3,kk+3)*A.*G(:,k+3,kk+3);
            end
            
            up=fft(up);
            % first order
            for l=1:M
                up(l)=exp(-1i*tau*betal(l))*up(l)+ delta_n1p(l)+ delta_n2p(l)+ delta_n3p(l);
            end
        end
        
        up=2*real(ifft(up));
        %% end
        us{fN,fe}=up;
        if fN>1
            err(fN-1,fe)=sqrt(h)*norm(us{fN,fe}-us{fN-1,fe},'fro');
        end
        %waitbar(fe*fN/14,ht);
    end%fN
end%fe
%close(ht);
t(fN,fe)=toc;
% x=linspace(-pi,pi-2*pi/M,M);
% plot(x,up)
%end



