%% 符号计算
%tips 1: 符号计算中可以嵌套共轭运算
%tips 2: 符号计算可以单独成文件先算，保存，载入直接赋值用symbolNPI_complex.m
function symbolNPI_complex(savefile,x)
%x=1->first-order NPI symbol integration
%x=2->first-order NPI symbol integration
%x=3->first-order NPI symbol integration
syms s w epsilo real;
%先判断x是否为整数
flagg=0;%用于指示是否已保存符号积分
if x==1 | x==2 | x==3
    %% first-order used funciton
    %p{k}(-s)=-conj(p{k}(s))
    p{1}=symfun(exp(-1i*s/epsilo^2)*s,[w,epsilo]);
    p{2}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*1i*w/epsilo^2),w,0,s),[w,epsilo]);
    p{3}=symfun(exp(-1i*s/epsilo^2)*int(exp(-2*1i*w/epsilo^2),w,0,s),[w,epsilo]);
    p{4}=symfun(exp(-1i*s/epsilo^2)*int(exp(4*1i*w/epsilo^2),w,0,s),[w,epsilo]);
    for k=1:4
        p{k}=subs(p{k},s,w);
    end
    if x==2 | x==3
        %% second-order used funciton
        %p{k}(-s)=conj(p{k}(s))
        p1{1}=symfun(exp(-1i*s/epsilo^2)*s^2/2,[w,epsilo]);
        p1{2}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*1i*w/epsilo^2)*w,w,0,s),[w,epsilo]);
        p1{3}=symfun(exp(-1i*s/epsilo^2)*int(exp(-2*1i*w/epsilo^2)*w,w,0,s),[w,epsilo]);
        p1{4}=symfun(exp(-1i*s/epsilo^2)*int(exp(4*1i*w/epsilo^2)*w,w,0,s),[w,epsilo]);
        for k=1:4
            p1{k}=subs(p1{k},s,w);
        end
        %1->plus,2->minus
        q{1,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(-conj(p{3})+p{4}),w,0,s),[w,epsilo]);
        q{1,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(-conj(p{4})+p{3}),w,0,s),[w,epsilo]);
        q{2,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(-conj(p{1})+p{2}),w,0,s),[w,epsilo]);
        q{2,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(-conj(p{2})+p{1}),w,0,s),[w,epsilo]);
        
        q{3,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(2*1i*w/epsilo^2)*(-conj(p{3})+p{4}),w,0,s),[w,epsilo]);
        q{3,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(-2*1i*w/epsilo^2)*(-conj(p{4})+p{3}),w,0,s),[w,epsilo]);
        q{4,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(2*1i*w/epsilo^2)*(-conj(p{1})+p{2}),w,0,s),[w,epsilo]);
        q{4,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(-2*1i*w/epsilo^2)*(-conj(p{2})+p{1}),w,0,s),[w,epsilo]);
        
        q{5,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(-2*1i*w/epsilo^2)*(-conj(p{3})+p{4}),w,0,s),[w,epsilo]);
        q{5,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(2*1i*w/epsilo^2)*(-conj(p{4})+p{3}),w,0,s),[w,epsilo]);
        q{6,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(-2*1i*w/epsilo^2)*(-conj(p{1})+p{2}),w,0,s),[w,epsilo]);
        q{6,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(2*1i*w/epsilo^2)*(-conj(p{2})+p{1}),w,0,s),[w,epsilo]);
        for j=1:6
            for k=1:2
                q{j,k}=subs(q{j,k},s,w);
            end
        end
        if x==3
            %%  third-order used function
            %p{k}(-s)=-conj(p{k}(s))
            %p_0,2/p_2,2/p_-2,2/p_4,2
            p2{1}=symfun(exp(-1i*s/epsilo^2)*s^3/3,[w,epsilo]);
            p2{2}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*1i*w/epsilo^2)*w^2,w,0,s),[w,epsilo]);
            p2{3}=symfun(exp(-1i*s/epsilo^2)*int(exp(-2*1i*w/epsilo^2)*w^2,w,0,s),[w,epsilo]);
            p2{4}=symfun(exp(-1i*s/epsilo^2)*int(exp(4*1i*w/epsilo^2)*w^2,w,0,s),[w,epsilo]);
            for k=1:4
                p2{k}=subs(p2{k},s,w);
            end
            %1->p,2->m
            %q_k,1,+/-, k=1,2,3,4,5,6
            q1{1,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(-conj(p{3})+p{4})*w,w,0,s),[w,epsilo]);
            q1{1,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(-conj(p{4})+p{3})*w,w,0,s),[w,epsilo]);
            q1{2,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(-conj(p{1})+p{2})*w,w,0,s),[w,epsilo]);
            q1{2,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(-conj(p{2})+p{1})*w,w,0,s),[w,epsilo]);
            
            q1{3,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(2*1i*w/epsilo^2)*(-conj(p{3})+p{4})*w,w,0,s),[w,epsilo]);
            q1{3,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(-2*1i*w/epsilo^2)*(-conj(p{4})+p{3})*w,w,0,s),[w,epsilo]);
            q1{4,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(2*1i*w/epsilo^2)*(-conj(p{1})+p{2})*w,w,0,s),[w,epsilo]);
            q1{4,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(-2*1i*w/epsilo^2)*(-conj(p{2})+p{1})*w,w,0,s),[w,epsilo]);
            
            q1{5,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(-2*1i*w/epsilo^2)*(-conj(p{3})+p{4})*w,w,0,s),[w,epsilo]);
            q1{5,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(2*1i*w/epsilo^2)*(-conj(p{4})+p{3})*w,w,0,s),[w,epsilo]);
            q1{6,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(-2*1i*w/epsilo^2)*(-conj(p{1})+p{2})*w,w,0,s),[w,epsilo]);
            q1{6,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(2*1i*w/epsilo^2)*(-conj(p{2})+p{1})*w,w,0,s),[w,epsilo]);
            %q_k,2,+/-, k=1,2,3,4,5,6
            q2{1,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(conj(p1{3})+p1{4}),w,0,s),[w,epsilo]);
            q2{1,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(conj(p1{4})+p1{3}),w,0,s),[w,epsilo]);
            q2{2,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(conj(p1{1})+p1{2}),w,0,s),[w,epsilo]);
            q2{2,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(conj(p1{2})+p1{1}),w,0,s),[w,epsilo]);
            
            q2{3,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(2*1i*w/epsilo^2)*(conj(p1{3})+p1{4}),w,0,s),[w,epsilo]);
            q2{3,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(-2*1i*w/epsilo^2)*(conj(p1{4})+p1{3}),w,0,s),[w,epsilo]);
            q2{4,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(2*1i*w/epsilo^2)*(conj(p1{1})+p1{2}),w,0,s),[w,epsilo]);
            q2{4,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(-2*1i*w/epsilo^2)*(conj(p1{2})+p1{1}),w,0,s),[w,epsilo]);
            
            q2{5,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(-2*1i*w/epsilo^2)*(conj(p1{3})+p1{4}),w,0,s),[w,epsilo]);
            q2{5,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(2*1i*w/epsilo^2)*(conj(p1{4})+p1{3}),w,0,s),[w,epsilo]);
            q2{6,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(-2*1i*w/epsilo^2)*(conj(p1{1})+p1{2}),w,0,s),[w,epsilo]);
            q2{6,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(2*1i*w/epsilo^2)*(conj(p1{2})+p1{1}),w,0,s),[w,epsilo]);
            %q_k,3,+/-, k=1,2,3,4,5,6
            q3{1,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(-conj(p{3})-p{4})*w,w,0,s),[w,epsilo]);
            q3{1,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(-conj(p{4})-p{3})*w,w,0,s),[w,epsilo]);
            q3{2,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(-conj(p{1})-p{2})*w,w,0,s),[w,epsilo]);
            q3{2,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(-conj(p{2})-p{1})*w,w,0,s),[w,epsilo]);
            
            q3{3,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(2*1i*w/epsilo^2)*(-conj(p{3})-p{4})*w,w,0,s),[w,epsilo]);
            q3{3,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(-2*1i*w/epsilo^2)*(-conj(p{4})-p{3})*w,w,0,s),[w,epsilo]);
            q3{4,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(2*1i*w/epsilo^2)*(-conj(p{1})-p{2})*w,w,0,s),[w,epsilo]);
            q3{4,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(-2*1i*w/epsilo^2)*(-conj(p{2})-p{1})*w,w,0,s),[w,epsilo]);
            
            q3{5,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(-2*1i*w/epsilo^2)*(-conj(p{3})-p{4})*w,w,0,s),[w,epsilo]);
            q3{5,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(2*1i*w/epsilo^2)*(-conj(p{4})-p{3})*w,w,0,s),[w,epsilo]);
            q3{6,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(-2*1i*w/epsilo^2)*(-conj(p{1})-p{2})*w,w,0,s),[w,epsilo]);
            q3{6,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(2*1i*w/epsilo^2)*(-conj(p{2})-p{1})*w,w,0,s),[w,epsilo]);
            
            
            q4{1,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(-conj(p1{3})+p1{4}),w,0,s),[w,epsilo]);
            q4{1,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(-conj(p1{4})+p1{3}),w,0,s),[w,epsilo]);
            q4{2,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(-conj(p1{1})+p1{2}),w,0,s),[w,epsilo]);
            q4{2,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(-conj(p1{2})+p1{1}),w,0,s),[w,epsilo]);
            
            q4{3,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(2*1i*w/epsilo^2)*(-conj(p1{3})+p1{4}),w,0,s),[w,epsilo]);
            q4{3,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(-2*1i*w/epsilo^2)*(-conj(p1{4})+p1{3}),w,0,s),[w,epsilo]);
            q4{4,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(2*1i*w/epsilo^2)*(-conj(p1{1})+p1{2}),w,0,s),[w,epsilo]);
            q4{4,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(-2*1i*w/epsilo^2)*(-conj(p1{2})+p1{1}),w,0,s),[w,epsilo]);
            
            q4{5,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(-2*1i*w/epsilo^2)*(-conj(p1{3})+p1{4}),w,0,s),[w,epsilo]);
            q4{5,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(2*1i*w/epsilo^2)*(-conj(p1{4})+p1{3}),w,0,s),[w,epsilo]);
            q4{6,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(-2*1i*w/epsilo^2)*(-conj(p1{1})+p1{2}),w,0,s),[w,epsilo]);
            q4{6,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(2*1i*w/epsilo^2)*(-conj(p1{2})+p1{1}),w,0,s),[w,epsilo]);
            for j=1:6
                for k=1:2
                    q1{j,k}=subs(q1{j,k},s,w);
                    q2{j,k}=subs(q2{j,k},s,w);
                    q3{j,k}=subs(q3{j,k},s,w);
                      q4{j,k}=subs(q4{j,k},s,w);
                end
            end
            
            %1->pp,2->pm,3->mp,4->mm
            %r_k,+/-,+/-, k=1,2,3,4,5,6,7,8 r的写法有所不同
            r{1,1}=symfun(exp(-1i*s/epsilo^2)*int((p{3}+-conj(p{4}))*(p{3}+-conj(p{4})),w,0,s),[w,epsilo]);
            r{1,2}=symfun(exp(-1i*s/epsilo^2)*int((p{3}+-conj(p{4}))*(-conj(p{3})+p{4}),w,0,s),[w,epsilo]);
            r{1,3}=symfun(exp(-1i*s/epsilo^2)*int((-conj(p{3})+p{4})*(p{3}+-conj(p{4})),w,0,s),[w,epsilo]);
            r{1,4}=symfun(exp(-1i*s/epsilo^2)*int((-conj(p{3})+p{4})*(-conj(p{3})+p{4}),w,0,s),[w,epsilo]);
            
            r{2,1}=symfun(exp(-1i*s/epsilo^2)*int((p{3}+-conj(p{4}))*(p{1}+-conj(p{2})),w,0,s),[w,epsilo]);
            r{2,2}=symfun(exp(-1i*s/epsilo^2)*int((p{3}+-conj(p{4}))*(-conj(p{1})+p{2}),w,0,s),[w,epsilo]);
            r{2,3}=symfun(exp(-1i*s/epsilo^2)*int((-conj(p{3})+p{4})*(p{1}+-conj(p{2})),w,0,s),[w,epsilo]);
            r{2,4}=symfun(exp(-1i*s/epsilo^2)*int((-conj(p{3})+p{4})*(-conj(p{1})+p{2}),w,0,s),[w,epsilo]);
            
            r{3,1}=symfun(exp(-1i*s/epsilo^2)*int((p{1}+-conj(p{2}))*(p{3}+-conj(p{4})),w,0,s),[w,epsilo]);
            r{3,2}=symfun(exp(-1i*s/epsilo^2)*int((p{1}+-conj(p{2}))*(-conj(p{3})+p{4}),w,0,s),[w,epsilo]);
            r{3,3}=symfun(exp(-1i*s/epsilo^2)*int((-conj(p{1})+p{2})*(p{3}+-conj(p{4})),w,0,s),[w,epsilo]);
            r{3,4}=symfun(exp(-1i*s/epsilo^2)*int((-conj(p{1})+p{2})*(-conj(p{3})+p{4}),w,0,s),[w,epsilo]);
            
            r{4,1}=symfun(exp(-1i*s/epsilo^2)*int((p{1}+-conj(p{2}))*(p{1}+-conj(p{2})),w,0,s),[w,epsilo]);
            r{4,2}=symfun(exp(-1i*s/epsilo^2)*int((p{1}+-conj(p{2}))*(-conj(p{1})+p{2}),w,0,s),[w,epsilo]);
            r{4,3}=symfun(exp(-1i*s/epsilo^2)*int((-conj(p{1})+p{2})*(p{1}+-conj(p{2})),w,0,s),[w,epsilo]);
            r{4,4}=symfun(exp(-1i*s/epsilo^2)*int((-conj(p{1})+p{2})*(-conj(p{1})+p{2}),w,0,s),[w,epsilo]);
            
            r{5,1}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*i*w/epsilo^2)*(p{3}+-conj(p{4}))*(p{3}+-conj(p{4})),w,0,s),[w,epsilo]);
            r{5,2}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*i*w/epsilo^2)*(p{3}+-conj(p{4}))*(-conj(p{3})+p{4}),w,0,s),[w,epsilo]);
            r{5,3}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*i*w/epsilo^2)*(-conj(p{3})+p{4})*(p{3}+-conj(p{4})),w,0,s),[w,epsilo]);
            r{5,4}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*i*w/epsilo^2)*(-conj(p{3})+p{4})*(-conj(p{3})+p{4}),w,0,s),[w,epsilo]);
            
            r{6,1}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*i*w/epsilo^2)*(p{3}+-conj(p{4}))*(p{1}+-conj(p{2})),w,0,s),[w,epsilo]);
            r{6,2}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*i*w/epsilo^2)*(p{3}+-conj(p{4}))*(-conj(p{1})+p{2}),w,0,s),[w,epsilo]);
            r{6,3}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*i*w/epsilo^2)*(-conj(p{3})+p{4})*(p{1}+-conj(p{2})),w,0,s),[w,epsilo]);
            r{6,4}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*i*w/epsilo^2)*(-conj(p{3})+p{4})*(-conj(p{1})+p{2}),w,0,s),[w,epsilo]);
            
            r{7,1}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*i*w/epsilo^2)*(p{1}+-conj(p{2}))*(p{3}+-conj(p{4})),w,0,s),[w,epsilo]);
            r{7,2}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*i*w/epsilo^2)*(p{1}+-conj(p{2}))*(-conj(p{3})+p{4}),w,0,s),[w,epsilo]);
            r{7,3}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*i*w/epsilo^2)*(-conj(p{1})+p{2})*(p{3}+-conj(p{4})),w,0,s),[w,epsilo]);
            r{7,4}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*i*w/epsilo^2)*(-conj(p{1})+p{2})*(-conj(p{3})+p{4}),w,0,s),[w,epsilo]);
            
            r{8,1}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*i*w/epsilo^2)*(p{1}+-conj(p{2}))*(p{1}+-conj(p{2})),w,0,s),[w,epsilo]);
            r{8,2}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*i*w/epsilo^2)*(p{1}+-conj(p{2}))*(-conj(p{1})+p{2}),w,0,s),[w,epsilo]);
            r{8,3}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*i*w/epsilo^2)*(-conj(p{1})+p{2})*(p{1}+-conj(p{2})),w,0,s),[w,epsilo]);
            r{8,4}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*i*w/epsilo^2)*(-conj(p{1})+p{2})*(-conj(p{1})+p{2}),w,0,s),[w,epsilo]);
            for j=1:8
                for k=1:4
                    r{j,k}=subs(r{j,k},s,w);
                end
            end
            %mk_j,+/-, k=1,2,3,j=1,2,3,4,5,6,7,8
            for k=1:6
                m1{k,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(q{k,1}+conj(q{k,2})),w,0,s),[w,epsilo]);
                m1{k,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(conj(q{k,1})+q{k,2}),w,0,s),[w,epsilo]);
                m2{k,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(2i*w/epsilo^2)*(q{k,1}+conj(q{k,2})),w,0,s),[w,epsilo]);
                m2{k,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(-2i*w/epsilo^2)*(conj(q{k,1})+q{k,2}),w,0,s),[w,epsilo]);
                m3{k,1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(-2i*w/epsilo^2)*(q{k,1}+conj(q{k,2})),w,0,s),[w,epsilo]);
                m3{k,2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(2i*w/epsilo^2)*(conj(q{k,1})+q{k,2}),w,0,s),[w,epsilo]);
            end
            for j=1:6
                for k=1:2
                    m1{j,k}=subs(m1{j,k},s,w);
                    m2{j,k}=subs(m2{j,k},s,w);
                    m3{j,k}=subs(m3{j,k},s,w);
                end
            end
            
             flagg=flagg+1;
        else
            disp('x=2');
        end%if
         flagg=flagg+1;
    else 
        disp('x=1');
    end
 flagg=flagg+1;
else
    disp('error');
end

if flagg==3
save(savefile,'p','p1','p2','q','q1','q2','q3','q4','r','m1','m2','m3');
elseif flagg ==2
save(savefile,'p','p1','q');
else 
save(savefile,'p');
end%function





