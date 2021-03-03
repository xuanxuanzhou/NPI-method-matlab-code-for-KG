%% 符号计算
%tips 1: 符号计算中可以嵌套共轭运算
%tips 2: 符号计算可以单独成文件先算，保存，载入直接赋值用symbolNPI0.m
syms s w epsilo real;
%first-order used funciton
%p{k}(-s)=-conj(p{k}(s))
p{1}=symfun(exp(-1i*s/epsilo^2)*s,[w,epsilo]);
p{2}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*1i*w/epsilo^2),w,0,s),[w,epsilo]);
p{3}=symfun(exp(-1i*s/epsilo^2)*int(exp(-2*1i*w/epsilo^2),w,0,s),[w,epsilo]);
p{4}=symfun(exp(-1i*s/epsilo^2)*int(exp(4*1i*w/epsilo^2),w,0,s),[w,epsilo]);
for k=1:4
    p{k}=subs(p{k},s,w);
end

%second-order used funciton
%p{k}(-s)=conj(p{k}(s))
p1{1}=symfun(exp(-1i*s/epsilo^2)*s^2/2,[w,epsilo]);
p1{2}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*1i*w/epsilo^2)*w,w,0,s),[w,epsilo]);
p1{3}=symfun(exp(-1i*s/epsilo^2)*int(exp(-2*1i*w/epsilo^2)*w,w,0,s),[w,epsilo]);
p1{4}=symfun(exp(-1i*s/epsilo^2)*int(exp(4*1i*w/epsilo^2)*w,w,0,s),[w,epsilo]);
for k=1:4
    p1{k}=subs(p1{k},s,w);
end
%1->p,2->m
q1{1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(-conj(p{3})+p{4}),w,0,s),[w,epsilo]);
q1{2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(-conj(p{4})+p{3}),w,0,s),[w,epsilo]);
q2{1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(-conj(p{1})+p{2}),w,0,s),[w,epsilo]);
q2{2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*(-conj(p{2})+p{1}),w,0,s),[w,epsilo]);

q3{1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(2*1i*w/epsilo^2)*(-conj(p{3})+p{4}),w,0,s),[w,epsilo]);
q3{2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(-2*1i*w/epsilo^2)*(-conj(p{4})+p{3}),w,0,s),[w,epsilo]);
q4{1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(2*1i*w/epsilo^2)*(-conj(p{1})+p{2}),w,0,s),[w,epsilo]);
q4{2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(-2*1i*w/epsilo^2)*(-conj(p{2})+p{1}),w,0,s),[w,epsilo]);

q5{1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(-2*1i*w/epsilo^2)*(-conj(p{3})+p{4}),w,0,s),[w,epsilo]);
q5{2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(2*1i*w/epsilo^2)*(-conj(p{4})+p{3}),w,0,s),[w,epsilo]);
q6{1}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(-2*1i*w/epsilo^2)*(-conj(p{1})+p{2}),w,0,s),[w,epsilo]);
q6{2}=symfun(int(exp(-1i*(s-w)/epsilo^2)*exp(2*1i*w/epsilo^2)*(-conj(p{2})+p{1}),w,0,s),[w,epsilo]);
for k=1:2
    q1{k}=subs(q1{k},s,w);
    q2{k}=subs(q2{k},s,w);
    q3{k}=subs(q3{k},s,w);
    q4{k}=subs(q4{k},s,w);
    q5{k}=subs(q5{k},s,w);
    q6{k}=subs(q6{k},s,w);
end

%third-order used function















