%% handbook
%the program is for the complex-value |u|^2u
%Klein-Gordon equation in the non-relativistic limit regime,
%first-order Nested Picard method
%with Fourier pseudospectral method,
%function: function [us,t]=newnpi1(flag)
%input:  the explanation of flag can be found below
%output: us is the numerical solution
%(the index of column is fe, the index of row is fN, t is the same as us),
%t is the cpu time
%%  function
function [us,t]=newnpi1(flag)
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
        M=8*16;                                                       
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
        for k=1:4                                                     
            %integral coefficient
            ps(k)=double(p{k}(tau,epsilon));
        end
        %computing
        %% computing
        for tk=1:N                                                    %notice the tk can not be used again below
            %delta^{n,1}_{+} delta_n1p
            ups=ifft(up);
            ums=ifft(um);
            F_1p=A.*fft(ups.^2.*conj(ums));
            F_1m=A.*fft(ums.^2.*conj(ups));
            F_2p=A.*fft((2*abs(ums).^2+abs(ups).^2).*ups);
            F_2m=A.*fft((2*abs(ups).^2+abs(ums).^2).*ums);
            
            delta_n1p=ps(1)*F_2p + ps(2)*F_2m+ ps(3)*F_1p + ps(4)*F_1m;
            delta_n1m=-conj(ps(1))*F_2m + -conj(ps(2))*F_2p+ -conj(ps(3))*F_1m + -conj(ps(4))*F_1p;
            
            ups=exp(-1i*tau*betal).*up+ delta_n1p;
            ums=exp(1i*tau*betal).*um+ delta_n1m;
            up=ups;
            um=ums;
        end
        %waitbar(fN*fe/21,ht);                                       %show the progress bar
        us{fN,fe}=ifft(up+um);                                              %save the numerical solution
    end%fN
    t(fN,fe)=toc;                                                         %save the cpu time
end%fe
%close(ht);                                                      %close the progress bar
end%function



