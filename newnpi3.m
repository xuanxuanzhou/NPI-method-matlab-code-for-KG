%% handbook
%the program is for the complex-value |u|^2u
%Klein-Gordon equation in the non-relativistic limit regime,
%third-order Nested Picard method
%with Fourier pseudospectral method,
%function: function [us,t]=newnpi3(flag)
%input: the explanation of flag can be found below
%output: us is the numerical solution
%(the index of column is fe, the index of row is fN, t is the same as us),
%t is the cpu time
%% function
%function [us,t]=newnpi3(flag)
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
        flag1=2;
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
            p2s(k)=double(p2{k}(tau,epsilon));
        end
        
        for j=1:6
            for k=1:2
                qs(j,k)=double(q{j,k}(tau,epsilon));
                q1s(j,k)=double(q1{j,k}(tau,epsilon));
                q2s(j,k)=double(q2{j,k}(tau,epsilon));
                q3s(j,k)=double(q3{j,k}(tau,epsilon));
                  q4s(j,k)=double(q4{j,k}(tau,epsilon));
            end
        end
        
        for j=1:8
            for k=1:4
                rs(j,k)=double(r{j,k}(tau,epsilon));
            end
        end
        
        for j=1:6
            for k=1:2
                m1s(j,k)=double(m1{j,k}(tau,epsilon));
                m2s(j,k)=double(m2{j,k}(tau,epsilon));
                m3s(j,k)=double(m3{j,k}(tau,epsilon));
            end
        end
        %% computing
        for tk=1:N                                                %notice the tk can not be used again below
            %delta^{n,1}_{+} delta_n1p
            ups=ifft(up);
            ums=ifft(um);
            ups1=ifft(B.*up);
            ums1=ifft(B.*um);
            ups2=ifft(B.^2.*up);
            ums2=ifft(B.^2.*um);
            %% first-order
            F_1p=A.*fft(ups.^2.*conj(ums));
            F_1m=A.*fft(ums.^2.*conj(ups));
            F_2p=A.*fft((2*abs(ums).^2+abs(ups).^2).*ups);
            F_2m=A.*fft((2*abs(ups).^2+abs(ums).^2).*ums);
            delta_n1p=ps(1)*F_2p + ps(2)*F_2m+ ps(3)*F_1p + ps(4)*F_1m;
            delta_n1m=-conj(ps(1))*F_2m + -conj(ps(2))*F_2p+ -conj(ps(3))*F_1m + -conj(ps(4))*F_1p;
             F_1p=ifft(F_1p);
             F_1m=ifft(F_1m);
             F_2p=ifft(F_2p);
             F_2m=ifft(F_2m);
            %% second-order
            G_1p=A.*B.*fft(ups.^2.*conj(ums));
            G_1m=A.*B.*fft(ums.^2.*conj(ups));
            G_2p=A.*B.*fft( (2*abs(ums).^2+ abs(ups).^2).*ups);
            G_2m=A.*B.*fft( (2*abs(ups).^2+ abs(ums).^2).*ums);
            
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

            delta_n2p=p1s(1)*(G_2p+G_4p) + p1s(2)*(G_2m+G_4m)+ p1s(3)*(G_1p+ G_3p) + p1s(4)*(G_1m+G_3m)...
                -tau*ps(1)*B.*fft(F_2p) - tau*ps(2)*B.*fft(F_2m)- tau*ps(3)*B.*fft(F_1p) - tau*ps(4)*B.*fft(F_1m)+...
                qs(1,1)*F_21p+ qs(1,2)*F_21m +qs(2,1)*F_22p+ qs(2,2)*F_22m +qs(3,1)*F_23p+ qs(3,2)*F_23m+...
                qs(4,1)*F_24p+ qs(4,2)*F_24m +qs(5,1)*F_25p+ qs(5,2)*F_25m +qs(6,1)*F_26p+ qs(6,2)*F_26m;
            
            delta_n2m=conj(p1s(1))*(G_2m-G_4m) + conj(p1s(2))*(G_2p- G_4p)+ conj(p1s(3))*(G_1m- G_3m) + conj(p1s(4))*(G_1p-G_3p)...
                -tau*conj(ps(1))*B.*fft(F_2m) - tau*conj(ps(2))*B.*fft(F_2p)- tau*conj(ps(3))*B.*fft(F_1m) - tau*conj(ps(4))*B.*fft(F_1p)+...
                conj(qs(1,1))*F_21m+ conj(qs(1,2))*F_21p + conj(qs(2,1))*F_22m+ conj(qs(2,2))*F_22p + conj(qs(3,1))*F_23m+ conj(qs(3,2))*F_23p+...
                conj(qs(4,1))*F_24m+ conj(qs(4,2))*F_24p + conj(qs(5,1))*F_25m+ conj(qs(5,2))*F_25p + conj(qs(6,1))*F_26m+ conj(qs(6,2))*F_26p;
            
             G_1p=ifft(G_1p);
             G_1m=ifft(G_1m);
             G_2p=ifft(G_2p);
             G_2m=ifft(G_2m);
             G_3p=ifft(G_3p);
             G_3m=ifft(G_3m);
             G_4p=ifft(G_4p);
             G_4m=ifft(G_4m);
            
            F_21p=ifft(F_21p);
            F_21m=ifft(F_21m);
            F_22p=ifft(F_22p);
            F_22m=ifft(F_22m);
            
            F_23p=ifft(F_23p);
            F_23m=ifft(F_23m);
            F_24p=ifft(F_24p);
            F_24m=ifft(F_24m);
            
            F_25p=ifft(F_25p);
            F_25m=ifft(F_25m);
            F_26p=ifft(F_26p);
            F_26m=ifft(F_26m);
            
            %% third-order(二阶向量，一阶向量均在物理空间)
            G_13p=A.*fft( 1/2*ups.^2.*conj(ums2)+ ups2.*ups.*conj(ums)  );
            G_13m=A.*fft( 1/2*ums.^2.*conj(ups2)+ ums2.*ums.*conj(ups)  );
            G_23p=A.*fft( ups.*conj(ums).*ums2 + (ups.*conj(ups)+ums.*conj(ums)).*ups2 +...
                1/2*ups.^2.*conj(ups2)+ conj(ums2).*ups.*ums   );
            G_23m=A.*fft( ums.*conj(ups).*ups2 + (ups.*conj(ups)+ums.*conj(ums)).*ums2 +...
                1/2*ums.^2.*conj(ums2)+ conj(ups2).*ums.*ups   );
            
            E_1p=A.*fft( -2*ups.*ups1.*conj(ums1) + conj(ums).*ups1.^2  );
            E_1m=A.*fft( -2*ums.*ums1.*conj(ups1) + conj(ups).*ums1.^2  );
            E_2p=A.*fft( 2*ups.*abs(ums1).^2 +2*ups.*abs(ups1).^2 -2*ums.*ups1.*conj(ums1) -2*conj(ums).*ups1.*ums1+ conj(ups).*ups1.^2  );
            E_2m=A.*fft( 2*ums.*abs(ups1).^2 +2*ums.*abs(ums1).^2 -2*ups.*ums1.*conj(ups1) -2*conj(ups).*ums1.*ups1+ conj(ums).*ums1.^2  );
            
            F_31p=A.*fft( 2*(ums.*conj(ums1)- ups.*conj(ups1) -conj(ups).*ups1+ conj(ums).* ums1 ).*F_1m- 2*(ups.*ums1- ums.*ups1).*conj(F_1p));
            F_31m=A.*fft( 2*(ums.*conj(ums1)- ups.*conj(ups1) -conj(ups).*ups1+ conj(ums).* ums1 ).*F_1p- 2*(ups.*ums1- ums.*ups1).*conj(F_1m));
            F_32p=A.*fft( 2*(ums.*conj(ums1)- ups.*conj(ups1) -conj(ups).*ups1+ conj(ums).* ums1 ).*F_2m- 2*(ups.*ums1- ums.*ups1).*conj(F_2p));
            F_32m=A.*fft( 2*(ums.*conj(ums1)- ups.*conj(ups1) -conj(ups).*ups1+ conj(ums).* ums1 ).*F_2p- 2*(ups.*ums1- ums.*ups1).*conj(F_2m));
            
            F_33p=A.*fft( 2*(-ums.*conj(ups1)+ conj(ups).*ums1 ).*F_1m-  2*(ums.*ums1).*conj(F_1p));
            F_33m=A.*fft( 2*(ups.*conj(ums1)- conj(ums).*ups1 ).*F_1p-  2*(-ups.*ups1).*conj(F_1m));
            F_34p=A.*fft( 2*(-ums.*conj(ups1)+ conj(ups).*ums1 ).*F_2m-  2*(ums.*ums1).*conj(F_2p));
            F_34m=A.*fft( 2*(ups.*conj(ums1)- conj(ums).*ups1 ).*F_2p-  2*(-ups.*ups1).*conj(F_2m));
            
            F_35p=A.*fft( 2*(ups.*conj(ums1)- conj(ums).*ups1 ).*F_1m-  2*(-ups.*ups1).*conj(F_1p));
            F_35m=A.*fft( 2*(-ums.*conj(ups1)+ conj(ups).*ums1 ).*F_1p-  2*(ums.*ums1).*conj(F_1m));
            F_36p=A.*fft( 2*(ups.*conj(ums1)- conj(ums).*ups1 ).*F_2m-  2*(-ups.*ups1).*conj(F_2p));
            F_36m=A.*fft( 2*(-ums.*conj(ups1)+ conj(ups).*ums1 ).*F_2p-  2*(ums.*ums1).*conj(F_2m));
            %%%%%%%%%%%%%%%%%%%%
            F_41p=A.*fft( 2*(abs(ups).^2+abs(ums).^2).*G_1m+ 2*ups.*ums.*conj(G_1p));
            F_41m=A.*fft( 2*(abs(ups).^2+abs(ums).^2).*G_1p+ 2*ups.*ums.*conj(G_1m));
            F_42p=A.*fft( 2*(abs(ups).^2+abs(ums).^2).*G_2m+ 2*ups.*ums.*conj(G_2p));
            F_42m=A.*fft( 2*(abs(ups).^2+abs(ums).^2).*G_2p+ 2*ups.*ums.*conj(G_2m));
            
            F_43p=A.*fft( 2*ums.*conj(ups).*G_1m+ ums.^2.*conj(G_1p));
            F_43m=A.*fft( 2*ups.*conj(ums).*G_1p+ ups.^2.*conj(G_1m));
            F_44p=A.*fft( 2*ums.*conj(ups).*G_2m+ ums.^2.*conj(G_2p));
            F_44m=A.*fft( 2*ups.*conj(ums).*G_2p+ ups.^2.*conj(G_2m));
            
            F_45p=A.*fft( 2*ups.*conj(ums).*G_1m+ ups.^2.*conj(G_1p));
            F_45m=A.*fft( 2*ums.*conj(ups).*G_1p+ ums.^2.*conj(G_1m));
            F_46p=A.*fft( 2*ups.*conj(ums).*G_2m+ ups.^2.*conj(G_2p));
            F_46m=A.*fft( 2*ums.*conj(ups).*G_2p+ ums.^2.*conj(G_2m));
            
            F_47p=A.*fft( 2*(abs(ups).^2+abs(ums).^2).*G_3m- 2*ups.*ums.*conj(G_3p));
            F_47m=A.*fft( 2*(abs(ups).^2+abs(ums).^2).*G_3p- 2*ups.*ums.*conj(G_3m));
            F_48p=A.*fft( 2*(abs(ups).^2+abs(ums).^2).*G_4m- 2*ups.*ums.*conj(G_4p));
            F_48m=A.*fft( 2*(abs(ups).^2+abs(ums).^2).*G_4p- 2*ups.*ums.*conj(G_4m));
            
            F_49p=A.*fft( 2*ums.*conj(ups).*G_3m- ums.^2.*conj(G_3p));
            F_49m=A.*fft( 2*ups.*conj(ums).*G_3p- ups.^2.*conj(G_3m));
            F_410p=A.*fft( 2*ums.*conj(ups).*G_4m- ums.^2.*conj(G_4p));
            F_410m=A.*fft( 2*ups.*conj(ums).*G_4p- ups.^2.*conj(G_4m));
            
            F_411p=A.*fft( 2*ups.*conj(ums).*G_3m- ups.^2.*conj(G_3p));
            F_411m=A.*fft( 2*ums.*conj(ups).*G_3p- ums.^2.*conj(G_3m));
            F_412p=A.*fft( 2*ups.*conj(ums).*G_4m- ups.^2.*conj(G_4p));
            F_412m=A.*fft( 2*ums.*conj(ups).*G_4p- ums.^2.*conj(G_4m));

            F_51p=A.*fft( 2*(abs(ups).^2+abs(ums).^2).*ifft(B.*fft(F_1m))+ 2*ups.*ums.*conj(ifft(B.*fft(F_1p))));
            F_51m=A.*fft( 2*(abs(ups).^2+abs(ums).^2).*ifft(B.*fft(F_1p))+2*ups.*ums.*conj(ifft(B.*fft(F_1m))));
            F_52p=A.*fft( 2*(abs(ups).^2+abs(ums).^2).*ifft(B.*fft(F_2m))+ 2*ups.*ums.*conj(ifft(B.*fft(F_2p))));
            F_52m=A.*fft( 2*(abs(ups).^2+abs(ums).^2).*ifft(B.*fft(F_2p))+ 2*ups.*ums.*conj(ifft(B.*fft(F_2m))));
            
            F_53p=A.*fft( 2*ums.*conj(ups).*ifft(B.*fft(F_1m))+ ums.^2.*conj(ifft(B.*fft(F_1p))));
            F_53m=A.*fft( 2*ups.*conj(ums).*ifft(B.*fft(F_1p))+ ups.^2.*conj(ifft(B.*fft(F_1m))));
            F_54p=A.*fft( 2*ums.*conj(ups).*ifft(B.*fft(F_2m))+ ums.^2.*conj(ifft(B.*fft(F_2p))));
            F_54m=A.*fft( 2*ups.*conj(ums).*ifft(B.*fft(F_2p))+ ups.^2.*conj(ifft(B.*fft(F_2m))));
            
            F_55p=A.*fft( 2*ups.*conj(ums).*ifft(B.*fft(F_1m))+ ups.^2.*conj(ifft(B.*fft(F_1p))));
            F_55m=A.*fft( 2*ums.*conj(ups).*ifft(B.*fft(F_1p))+ ums.^2.*conj(ifft(B.*fft(F_1m))));
            F_56p=A.*fft( 2*ups.*conj(ums).*ifft(B.*fft(F_2m))+ ups.^2.*conj(ifft(B.*fft(F_2p))));
            F_56m=A.*fft( 2*ums.*conj(ups).*ifft(B.*fft(F_2p))+ ums.^2.*conj(ifft(B.*fft(F_2m))));
            
            F_61_1p=A.*fft( 2*(abs(ups).^2+abs(ums).^2).*F_21p+ 2*ups.*ums.*conj(F_21m));
            F_61_1m=A.*fft( 2*(abs(ups).^2+abs(ums).^2).*F_21m+ 2*ups.*ums.*conj(F_21p));
            F_61_2p=A.*fft( 2*(abs(ups).^2+abs(ums).^2).*F_22p+ 2*ups.*ums.*conj(F_22m));
            F_61_2m=A.*fft( 2*(abs(ups).^2+abs(ums).^2).*F_22m+ 2*ups.*ums.*conj(F_22p));
            F_61_3p=A.*fft( 2*(abs(ups).^2+abs(ums).^2).*F_23p+ 2*ups.*ums.*conj(F_23m));
            F_61_3m=A.*fft( 2*(abs(ups).^2+abs(ums).^2).*F_23m+ 2*ups.*ums.*conj(F_23p));
            F_61_4p=A.*fft( 2*(abs(ups).^2+abs(ums).^2).*F_24p+ 2*ups.*ums.*conj(F_24m));
            F_61_4m=A.*fft( 2*(abs(ups).^2+abs(ums).^2).*F_24m+ 2*ups.*ums.*conj(F_24p));
            F_61_5p=A.*fft( 2*(abs(ups).^2+abs(ums).^2).*F_25p+ 2*ups.*ums.*conj(F_25m));
            F_61_5m=A.*fft( 2*(abs(ups).^2+abs(ums).^2).*F_25m+ 2*ups.*ums.*conj(F_25p));
            F_61_6p=A.*fft( 2*(abs(ups).^2+abs(ums).^2).*F_26p+ 2*ups.*ums.*conj(F_26m));
            F_61_6m=A.*fft( 2*(abs(ups).^2+abs(ums).^2).*F_26m+ 2*ups.*ums.*conj(F_26p));
            
            F_62_1p=A.*fft( 2*ums.*conj(ups).*F_21p+ ums.^2.*conj(F_21m));
            F_62_1m=A.*fft( 2*ups.*conj(ums).*F_21m+ ups.^2.*conj(F_21p));
            F_62_2p=A.*fft( 2*ums.*conj(ups).*F_22p+ ums.^2.*conj(F_22m));
            F_62_2m=A.*fft( 2*ups.*conj(ums).*F_22m+ ups.^2.*conj(F_22p));
            F_62_3p=A.*fft( 2*ums.*conj(ups).*F_23p+ ums.^2.*conj(F_23m));
            F_62_3m=A.*fft( 2*ups.*conj(ums).*F_23m+ ups.^2.*conj(F_23p));
            F_62_4p=A.*fft( 2*ums.*conj(ups).*F_24p+ ums.^2.*conj(F_24m));
            F_62_4m=A.*fft( 2*ups.*conj(ums).*F_24m+ ups.^2.*conj(F_24p));
            F_62_5p=A.*fft( 2*ums.*conj(ups).*F_25p+ ums.^2.*conj(F_25m));
            F_62_5m=A.*fft( 2*ups.*conj(ums).*F_25m+ ups.^2.*conj(F_25p));
            F_62_6p=A.*fft( 2*ums.*conj(ups).*F_26p+ ums.^2.*conj(F_26m));
            F_62_6m=A.*fft( 2*ups.*conj(ums).*F_26m+ ups.^2.*conj(F_26p));
            
            F_63_1p=A.*fft( 2*ups.*conj(ums).*F_21p+ ups.^2.*conj(F_21m));
            F_63_1m=A.*fft( 2*ums.*conj(ups).*F_21m+ ums.^2.*conj(F_21p));
            F_63_2p=A.*fft( 2*ups.*conj(ums).*F_22p+ ups.^2.*conj(F_22m));
            F_63_2m=A.*fft( 2*ums.*conj(ups).*F_22m+ ums.^2.*conj(F_22p));
            F_63_3p=A.*fft( 2*ups.*conj(ums).*F_23p+ ups.^2.*conj(F_23m));
            F_63_3m=A.*fft( 2*ums.*conj(ups).*F_23m+ ums.^2.*conj(F_23p));
            F_63_4p=A.*fft( 2*ups.*conj(ums).*F_24p+ ups.^2.*conj(F_24m));
            F_63_4m=A.*fft( 2*ums.*conj(ups).*F_24m+ ums.^2.*conj(F_24p));
            F_63_5p=A.*fft( 2*ups.*conj(ums).*F_25p+ ups.^2.*conj(F_25m));
            F_63_5m=A.*fft( 2*ums.*conj(ups).*F_25m+ ums.^2.*conj(F_25p));
            F_63_6p=A.*fft( 2*ups.*conj(ums).*F_26m+ ups.^2.*conj(F_26m));
            F_63_6m=A.*fft( 2*ums.*conj(ups).*F_26m+ ums.^2.*conj(F_26p));
            
            H_1ppp=A.*fft( -2*ups.*F_1p.*conj(F_1m) + conj(ums).*F_1p.*F_1p);
            H_1ppm=A.*fft( -2*ups.*F_1p.*conj(F_1p) + conj(ums).*F_1p.*F_1m);
            H_1pmp=A.*fft( -2*ups.*F_1m.*conj(F_1m) + conj(ums).*F_1m.*F_1p);
            H_1pmm=A.*fft( -2*ups.*F_1m.*conj(F_1p) + conj(ums).*F_1m.*F_1m);
            H_1mpp=A.*fft( -2*ums.*F_1p.*conj(F_1m) + conj(ups).*F_1p.*F_1p);
            H_1mpm=A.*fft( -2*ums.*F_1p.*conj(F_1p) + conj(ups).*F_1p.*F_1m);
            H_1mmp=A.*fft( -2*ums.*F_1m.*conj(F_1m) + conj(ups).*F_1m.*F_1p);
            H_1mmm=A.*fft( -2*ums.*F_1m.*conj(F_1p) + conj(ups).*F_1m.*F_1m);
            
            H_2ppp=A.*fft( -2*ups.*F_1p.*conj(F_2m) + conj(ums).*F_1p.*F_2p);
            H_2ppm=A.*fft( -2*ups.*F_1p.*conj(F_2p) + conj(ums).*F_1p.*F_2m);
            H_2pmp=A.*fft( -2*ups.*F_1m.*conj(F_2m) + conj(ums).*F_1m.*F_2p);
            H_2pmm=A.*fft( -2*ups.*F_1m.*conj(F_2p) + conj(ums).*F_1m.*F_2m);
            H_2mpp=A.*fft( -2*ums.*F_1p.*conj(F_2m) + conj(ups).*F_1p.*F_2p);
            H_2mpm=A.*fft( -2*ums.*F_1p.*conj(F_2p) + conj(ups).*F_1p.*F_2m);
            H_2mmp=A.*fft( -2*ums.*F_1m.*conj(F_2m) + conj(ups).*F_1m.*F_2p);
            H_2mmm=A.*fft( -2*ums.*F_1m.*conj(F_2p) + conj(ups).*F_1m.*F_2m);
            
            H_3ppp=A.*fft( -2*ups.*F_2p.*conj(F_1m) + conj(ums).*F_2p.*F_1p);
            H_3ppm=A.*fft( -2*ups.*F_2p.*conj(F_1p) + conj(ums).*F_2p.*F_1m);
            H_3pmp=A.*fft( -2*ups.*F_2m.*conj(F_1m) + conj(ums).*F_2m.*F_1p);
            H_3pmm=A.*fft( -2*ups.*F_2m.*conj(F_1p) + conj(ums).*F_2m.*F_1m);
            H_3mpp=A.*fft( -2*ums.*F_2p.*conj(F_1m) + conj(ups).*F_2p.*F_1p);
            H_3mpm=A.*fft( -2*ums.*F_2p.*conj(F_1p) + conj(ups).*F_2p.*F_1m);
            H_3mmp=A.*fft( -2*ums.*F_2m.*conj(F_1m) + conj(ups).*F_2m.*F_1p);
            H_3mmm=A.*fft( -2*ums.*F_2m.*conj(F_1p) + conj(ups).*F_2m.*F_1m);
            
            H_4ppp=A.*fft( -2*ups.*F_2p.*conj(F_2m) + conj(ums).*F_2p.*F_2p);
            H_4ppm=A.*fft( -2*ups.*F_2p.*conj(F_2p) + conj(ums).*F_2p.*F_2m);
            H_4pmp=A.*fft( -2*ups.*F_2m.*conj(F_2m) + conj(ums).*F_2m.*F_2p);
            H_4pmm=A.*fft( -2*ups.*F_2m.*conj(F_2p) + conj(ums).*F_2m.*F_2m);
            H_4mpp=A.*fft( -2*ums.*F_2p.*conj(F_2m) + conj(ups).*F_2p.*F_2p);
            H_4mpm=A.*fft( -2*ums.*F_2p.*conj(F_2p) + conj(ups).*F_2p.*F_2m);
            H_4mmp=A.*fft( -2*ums.*F_2m.*conj(F_2m) + conj(ups).*F_2m.*F_2p);
            H_4mmm=A.*fft( -2*ums.*F_2m.*conj(F_2p) + conj(ups).*F_2m.*F_2m);

            delta_n3p=tau^2*1/2*B.^2.*(   ps(1)*fft(F_2p) + ps(2)*fft(F_2m)+ ps(3)*fft(F_1p) + ps(4)*fft(F_1m)  )...
                +1/2*B.^2.*( p2s(1)*fft(F_2p) + p2s(2)*fft(F_2m)+ p2s(3)*fft(F_1p) + p2s(4)*fft(F_1m)  )...
                -tau*B.^2.*( p1s(1)*fft(F_2p) + p1s(2)*fft(F_2m)+ p1s(3)*fft(F_1p) + p1s(4)*fft(F_1m)  )...
                -tau*B.*( p1s(1)*fft(G_4p) + p1s(2)*fft(G_4m)+ p1s(3)*fft(+ G_3p) + p1s(4)*fft(G_3m)   )...
                +B.*( p2s(1)*fft(G_4p) + p2s(2)*fft(G_4m)+ p2s(3)*fft(+ G_3p) + p2s(4)*fft(G_3m)   )...
                  +( p2s(1)*E_2p + p2s(2)*E_2m+ p2s(3)*E_1p + p2s(4)*E_1m  )...
                   + p2s(1)*G_23p + p2s(2)*G_23m+ p2s(3)*G_13p + p2s(4)*G_13m...
                 -tau*B.*(  qs(1,1)*fft(F_21p)+ qs(1,2)*fft(F_21m) +qs(2,1)*fft(F_22p)+ qs(2,2)*fft(F_22m) +qs(3,1)*fft(F_23p)+ qs(3,2)*fft(F_23m)+...
                qs(4,1)*fft(F_24p)+ qs(4,2)*fft(F_24m) +qs(5,1)*fft(F_25p)+ qs(5,2)*fft(F_25m) +qs(6,1)*fft(F_26p)+ qs(6,2)*fft(F_26m)  )...
                   +B.* ( q1s(1,1)*fft(F_21p)+ q1s(1,2)*fft(F_21m) +q1s(2,1)*fft(F_22p)+ q1s(2,2)*fft(F_22m) +q1s(3,1)*fft(F_23p)+ q1s(3,2)*fft(F_23m)+...
                q1s(4,1)*fft(F_24p)+ q1s(4,2)*fft(F_24m) +q1s(5,1)*fft(F_25p)+ q1s(5,2)*fft(F_25m) +q1s(6,1)*fft(F_26p)+ q1s(6,2)*fft(F_26m))...
               +(  q1s(1,1)*F_31p+ q1s(1,2)*F_31m +q1s(2,1)*F_32p+ q1s(2,2)*F_32m +q1s(3,1)*F_33p+ q1s(3,2)*F_33m+...
                    q1s(4,1)*F_34p+ q1s(4,2)*F_34m +q1s(5,1)*F_35p+ q1s(5,2)*F_35m +q1s(6,1)*F_36p+ q1s(6,2)*F_36m  )...
                +rs(1,1)*H_1ppp +rs(1,2)*H_1ppm+  rs(1,3)*H_1pmp+ rs(1,4)*H_1pmm+ rs(2,1)*H_2ppp +rs(2,2)*H_2ppm+  rs(2,3)*H_2pmp+ rs(2,4)*H_2pmm...
                +rs(3,1)*H_3ppp +rs(3,2)*H_3ppm+  rs(3,3)*H_3pmp+ rs(3,4)*H_3pmm+ rs(4,1)*H_4ppp +rs(4,2)*H_4ppm+  rs(4,3)*H_4pmp+ rs(4,4)*H_4pmm...
                +rs(5,1)*H_1mpp+rs(5,2)*H_1mpm+ rs(5,3)*H_1mmp+ rs(5,4)*H_1mmm+ rs(6,1)*H_2mpp +rs(6,2)*H_2mpm+  rs(6,3)*H_2mmp+ rs(6,4)*H_2mmm...
                +rs(7,1)*H_3mpp+rs(7,2)*H_3mpm+ rs(7,3)*H_3mmp+ rs(7,4)*H_3mmm+ rs(8,1)*H_4mpp +rs(8,2)*H_4mpm+  rs(8,3)*H_4mmp+ rs(8,4)*H_4mmm...
                +m1s(1,1)*F_61_1p+  m1s(2,1)*F_61_2p+ m1s(3,1)*F_61_3p+ m1s(4,1)*F_61_4p+ m1s(5,1)*F_61_5p+ m1s(6,1)*F_61_6p+...
                +m2s(1,1)*F_62_1p+  m2s(2,1)*F_62_2p+ m2s(3,1)*F_62_3p+ m2s(4,1)*F_62_4p+ m2s(5,1)*F_62_5p+ m2s(6,1)*F_62_6p+...
                +m3s(1,1)*F_63_1p+  m3s(2,1)*F_63_2p+ m3s(3,1)*F_63_3p+ m3s(4,1)*F_63_4p+ m3s(5,1)*F_63_5p+ m3s(6,1)*F_63_6p+...
                +m1s(1,2)*F_61_1m+  m1s(2,2)*F_61_2m+ m1s(3,2)*F_61_3m+ m1s(4,2)*F_61_4m+ m1s(5,2)*F_61_5m+ m1s(6,2)*F_61_6m+...
                +m2s(1,2)*F_62_1m+  m2s(2,2)*F_62_2m+ m2s(3,2)*F_62_3m+ m2s(4,2)*F_62_4m+ m2s(5,2)*F_62_5m+ m2s(6,2)*F_62_6m+...
                +m3s(1,2)*F_63_1m+  m3s(2,2)*F_63_2m+ m3s(3,2)*F_63_3m+ m3s(4,2)*F_63_4m+ m3s(5,2)*F_63_5m+ m3s(6,2)*F_63_6m...
                  +q2s(1,1)*F_41p+ q2s(1,2)*F_41m +q2s(2,1)*F_42p+ q2s(2,2)*F_42m +q2s(3,1)*F_43p+ q2s(3,2)*F_43m+...
                  q2s(4,1)*F_44p+ q2s(4,2)*F_44m +q2s(5,1)*F_45p+ q2s(5,2)*F_45m +q2s(6,1)*F_46p+ q2s(6,2)*F_46m...
                +q3s(1,1)*F_51p+ q3s(1,2)*F_51m +q3s(2,1)*F_52p+ q3s(2,2)*F_52m +q3s(3,1)*F_53p+ q3s(3,2)*F_53m+...
                  q3s(4,1)*F_54p+ q3s(4,2)*F_54m +q3s(5,1)*F_55p+q3s(5,2)*F_55m +q3s(6,1)*F_56p+ q3s(6,2)*F_56m...
                +q4s(1,1)*F_47p+ q4s(1,2)*F_47m +q4s(2,1)*F_48p+ q4s(2,2)*F_48m +q4s(3,1)*F_49p+ q4s(3,2)*F_49m+...
                  +q4s(4,1)*F_410p+ q4s(4,2)*F_410m +q4s(5,1)*F_411p+ q4s(5,2)*F_411m +q4s(6,1)*F_412p+ q4s(6,2)*F_412m;
              
           delta_n3m=tau^2*1/2*B.^2.*(  -conj(ps(1))*fft(F_2m) + -conj(ps(2))*fft(F_2p)+ -conj(ps(3))*fft(F_1m) + -conj(ps(4))*fft(F_1p)  )+...
               -1/2*B.^2.*( conj(p2s(1))*fft(F_2m) + conj(p2s(2))*fft(F_2p)+ conj(p2s(3))*fft(F_1m) + conj(p2s(4))*fft(F_1p)  )...
                +tau*B.^2.*( conj(p1s(1))*fft(F_2m) + conj(p1s(2))*fft(F_2p)+ conj(p1s(3))*fft(F_1m) + conj(p1s(4))*fft(F_1p)  )...
                -tau*B.*( conj(p1s(1))*fft(G_4m) + conj(p1s(2))*fft(G_4p)+ conj(p1s(3))*fft(+ G_3m) + conj(p1s(4))*fft(G_3p)   )...
                +B.*( conj(p2s(1))*fft(G_4m) + conj(p2s(2))*fft(G_4p)+ conj(p2s(3))*fft(+ G_3m) + conj(p2s(4))*fft(G_3p)   )...
                -( conj(p2s(1))*E_2m + conj(p2s(2))*E_2p+ conj(p2s(3))*E_1m + conj(p2s(4))*E_1p  )...
              -( conj(p2s(1))*G_23m + conj(p2s(2))*G_23p+ conj(p2s(3))*G_13m + conj(p2s(4))*G_13p   )...
                +tau*B.* ( conj(qs(1,1))*fft(F_21m)+ conj(qs(1,2))*fft(F_21p) + conj(qs(2,1))*fft(F_22m)+ conj(qs(2,2))*fft(F_22p) + conj(qs(3,1))*fft(F_23m)+ conj(qs(3,2))*fft(F_23p)+...
                conj(qs(4,1))*fft(F_24m)+ conj(qs(4,2))*fft(F_24p) + conj(qs(5,1))*fft(F_25m)+ conj(qs(5,2))*fft(F_25p) +conj(qs(6,1))*fft(F_26m)+ conj(qs(6,2))*fft(F_26p) )...
             -B.*( conj(q1s(1,1))*fft(F_21m)+ conj(q1s(1,2))*fft(F_21p) + conj(q1s(2,1))*fft(F_22m)+ conj(q1s(2,2))*fft(F_22p) + conj(q1s(3,1))*fft(F_23m)+ conj(q1s(3,2))*fft(F_23p)+...
                conj(q1s(4,1))*fft(F_24m)+ conj(q1s(4,2))*fft(F_24p) + conj(q1s(5,1))*fft(F_25m)+ conj(q1s(5,2))*fft(F_25p) + conj(q1s(6,1))*fft(F_26m)+ conj(q1s(6,2))*fft(F_26p) )...
                  +(  conj(q1s(1,1))*F_31m+ conj(q1s(1,2))*F_31p + conj(q1s(2,1))*F_32m+ conj(q1s(2,2))*F_32p + conj(q1s(3,1))*F_33m+ conj(q1s(3,2))*F_33p+...
                    conj(q1s(4,1))*F_34m+ conj(q1s(4,2))*F_34p + conj(q1s(5,1))*F_35m+ conj(q1s(5,2))*F_35p + conj(q1s(6,1))*F_36m+ conj(q1s(6,2))*F_36p  )...
             -( conj(rs(1,1))*H_1mpp + conj(rs(1,2))*H_1mpm+  conj(rs(1,3))*H_1mmp+ conj(rs(1,4))*H_1mmm+ conj(rs(2,1))*H_2mpp + conj(rs(2,2))*H_2mpm+  conj(rs(2,3))*H_2mmp+ conj(rs(2,4))*H_2mmm...
                +conj(rs(3,1))*H_3mpp +conj(rs(3,2))*H_3mpm+  conj(rs(3,3))*H_3mmp+ conj(rs(3,4))*H_3mmm+ conj(rs(4,1))*H_4mpp +conj(rs(4,2))*H_4mpm+ conj( rs(4,3))*H_4mmp+ conj(rs(4,4))*H_4mmm...
                +conj(rs(5,1))*H_1ppp+conj(rs(5,2))*H_1ppm+ conj(rs(5,3))*H_1pmp+ conj(rs(5,4))*H_1pmm+ conj(rs(6,1))*H_2ppp +conj(rs(6,2))*H_2ppm+  conj(rs(6,3))*H_2pmp+ conj(rs(6,4))*H_2pmm...
                +conj(rs(7,1))*H_3ppp+conj(rs(7,2))*H_3ppm+ conj(rs(7,3))*H_3pmp+ conj(rs(7,4))*H_3pmm+ conj(rs(8,1))*H_4ppp +conj(rs(8,2))*H_4ppm+  conj(rs(8,3))*H_4pmp+ conj(rs(8,4))*H_4pmm  )...
                 -(+conj(q2s(1,1))*F_41m+ conj(q2s(1,2))*F_41p + conj(q2s(2,1))*F_42m+ conj(q2s(2,2))*F_42p + conj(q2s(3,1))*F_43m+ conj(q2s(3,2))*F_43p+...
                  conj(q2s(4,1))*F_44m+ conj(q2s(4,2))*F_44p + conj(q2s(5,1))*F_45m+ conj(q2s(5,2))*F_45p + conj(q2s(6,1))*F_46m+ conj(q2s(6,2))*F_46p )...
                  -( conj(q4s(1,1))*F_47m+ conj(q4s(1,2))*F_47p +conj(q4s(2,1))*F_48m+ conj(q4s(2,2))*F_48p +conj(q4s(3,1))*F_49m+ conj(q4s(3,2))*F_49p+...
                  conj(q4s(4,1))*F_410m+ conj(q4s(4,2))*F_410p +conj(q4s(5,1))*F_411m+ conj(q4s(5,2))*F_411p +conj(q4s(6,1))*F_412m+ conj(q4s(6,2))*F_412p  )...
                  -( conj(q3s(1,1))*F_51m+ conj(q3s(1,2))*F_51p + conj(q3s(2,1))*F_52m+ conj(q3s(2,2))*F_52p + conj(q3s(3,1))*F_53m+ conj(q3s(3,2))*F_53p+...
                  conj(q3s(4,1))*F_54m+ conj(q3s(4,2))*F_54p + conj(q3s(5,1))*F_55m+ conj(q3s(5,2))*F_55p + conj(q3s(6,1))*F_56m+ conj(q3s(6,2))*F_56p )...
                -(  conj(m1s(1,1))*F_61_1m+  conj(m1s(2,1))*F_61_2m+ conj(m1s(3,1))*F_61_3m+ conj(m1s(4,1))*F_61_4m+ conj(m1s(5,1))*F_61_5m+ conj(m1s(6,1))*F_61_6m+...
                +conj(m2s(1,1))*F_62_1m+  conj(m2s(2,1))*F_62_2m+ conj(m2s(3,1))*F_62_3m+ conj(m2s(4,1))*F_62_4m+ conj(m2s(5,1))*F_62_5m+ conj(m2s(6,1))*F_62_6m+...
                +conj(m3s(1,1))*F_63_1m+  conj(m3s(2,1))*F_63_2m+ conj(m3s(3,1))*F_63_3m+ conj(m3s(4,1))*F_63_4m+ conj(m3s(5,1))*F_63_5m+ conj(m3s(6,1))*F_63_6m+...
                +conj(m1s(1,2))*F_61_1p+  conj(m1s(2,2))*F_61_2p+ conj(m1s(3,2))*F_61_3p+ conj(m1s(4,2))*F_61_4p+ conj(m1s(5,2))*F_61_5p+ conj(m1s(6,2))*F_61_6p+...
                +conj(m2s(1,2))*F_62_1p+  conj(m2s(2,2))*F_62_2p+ conj(m2s(3,2))*F_62_3p+ conj(m2s(4,2))*F_62_4p+ conj(m2s(5,2))*F_62_5p+ conj(m2s(6,2))*F_62_6p+...
                +conj(m3s(1,2))*F_63_1p+  conj(m3s(2,2))*F_63_2p+ conj(m3s(3,2))*F_63_3p+ conj(m3s(4,2))*F_63_4p+ conj(m3s(5,2))*F_63_5p+ conj(m3s(6,2))*F_63_6p   );
            
            ups=exp(-1i*tau*betal).*up+ delta_n1p+ delta_n2p+ delta_n3p;
            ums=exp(1i*tau*betal).*um+ delta_n1m+ delta_n2m+ delta_n3m;
            up=ups;
            um=ums;
        end
        % waitbar(fN*fe/21,ht);                                       %show the progress bar
        us{fN,fe}=ifft(up+um);
        if fN>1
            err(fN-1,fe)=sqrt(h)*norm(us{fN,fe}-us{fN-1,fe},'fro');
        end
    end%fN
    t(fN,fe)=toc;                                                         %save the cpu time
end%fe
%end%function