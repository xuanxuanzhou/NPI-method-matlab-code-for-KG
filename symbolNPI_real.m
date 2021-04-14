%% 符号计算
%tips 1: 符号计算中可以嵌套共轭运算
%tips 2: 符号计算可以单独成文件先算，保存，载入直接赋值用symbolNPI.m
syms s w epsilo real;
%first-order used funciton
p{1}=symfun(exp(-1i*s/epsilo^2)*int(exp(3*1i*w/epsilo^2),w,0,s),[w,epsilo]);
p{2}=symfun(exp(-1i*s/epsilo^2)*int(exp(1*1i*w/epsilo^2),w,0,s),[w,epsilo]);
p{3}=symfun(exp(-1i*s/epsilo^2)*int(exp(-1*1i*w/epsilo^2),w,0,s),[w,epsilo]);
for k=1:3
    p{k}=subs(p{k},s,w);
end
%second-order used funciton
p1{1}=symfun(exp(-1i*s/epsilo^2)*int(exp(3*1i*w/epsilo^2)*w,w,0,s),[w,epsilo]);
p1{2}=symfun(exp(-1i*s/epsilo^2)*int(exp(1*1i*w/epsilo^2)*w,w,0,s),[w,epsilo]);
p1{3}=symfun(exp(-1i*s/epsilo^2)*int(exp(-1*1i*w/epsilo^2)*w,w,0,s),[w,epsilo]);
for k=1:3
    p1{k}=subs(p1{k},s,w);
end
for k=1:3
    q2{k}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*1i*w/epsilo^2)*p{k},w,0,s),[w,epsilo]);
    q2c{k}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*1i*w/epsilo^2)*conj(p{k}),w,0,s),[w,epsilo]);
    q0{k}=symfun(exp(-1i*s/epsilo^2)*int(exp(0*1i*w/epsilo^2)*p{k},w,0,s),[w,epsilo]);
    q0c{k}=symfun(exp(-1i*s/epsilo^2)*int(exp(0*1i*w/epsilo^2)*conj(p{k}),w,0,s),[w,epsilo]);

end
for k=1:3
    q2{k}=subs(q2{k},s,w);
    q2c{k}=subs(q2c{k},s,w);
    q0{k}=subs(q0{k},s,w);
    q0c{k}=subs(q0c{k},s,w);
end
%third-order used function
p2{1}=symfun(exp(-1i*s/epsilo^2)*int(exp(3*1i*w/epsilo^2)*w*w,w,0,s),[w,epsilo]);
p2{2}=symfun(exp(-1i*s/epsilo^2)*int(exp(1*1i*w/epsilo^2)*w*w,w,0,s),[w,epsilo]);
p2{3}=symfun(exp(-1i*s/epsilo^2)*int(exp(-1*1i*w/epsilo^2)*w*w,w,0,s),[w,epsilo]);
for k=1:3
    p2{k}=subs(p2{k},s,w);
end
for k=1:3
    q21{k}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*1i*w/epsilo^2)*p{k}*w,w,0,s),[w,epsilo]);
    q2c1{k}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*1i*w/epsilo^2)*conj(p{k})*w,w,0,s),[w,epsilo]);
    q01{k}=symfun(exp(-1i*s/epsilo^2)*int(exp(0*1i*w/epsilo^2)*p{k}*w,w,0,s),[w,epsilo]);
    q0c1{k}=symfun(exp(-1i*s/epsilo^2)*int(exp(0*1i*w/epsilo^2)*conj(p{k})*w,w,0,s),[w,epsilo]);
end
for k=1:3
    q21{k}=subs(q21{k},s,w);
    q2c1{k}=subs(q2c1{k},s,w);
    q01{k}=subs(q01{k},s,w);
    q0c1{k}=subs(q0c1{k},s,w);
end
for k=1:3
    for kk=1:3
         r{k,kk}=symfun(exp(-1i*s/epsilo^2)*int(exp(1i*w/epsilo^2)*p{k}*p{kk},w,0,s),[w,epsilo]);
          r{k+3,kk}=symfun(exp(-1i*s/epsilo^2)*int(exp(1i*w/epsilo^2)*conj(p{k})*p{kk},w,0,s),[w,epsilo]);
           r{k,kk+3}=symfun(exp(-1i*s/epsilo^2)*int(exp(1i*w/epsilo^2)*p{k}*conj(p{kk}),w,0,s),[w,epsilo]);
            r{k+3,kk+3}=symfun(exp(-1i*s/epsilo^2)*int(exp(1i*w/epsilo^2)*conj(p{k}*p{kk}),w,0,s),[w,epsilo]);
    end
end
for k=1:6
    for kk=1:6
         r{k,kk}=subs(r{k,kk},s,w);
    end
end

for k=1:3
    m2{k}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*1i*w/epsilo^2)*p1{k},w,0,s),[w,epsilo]);
    m2c{k}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*1i*w/epsilo^2)*conj(p1{k}),w,0,s),[w,epsilo]);
    m0{k}=symfun(exp(-1i*s/epsilo^2)*int(exp(0*1i*w/epsilo^2)*p1{k},w,0,s),[w,epsilo]);
    m0c{k}=symfun(exp(-1i*s/epsilo^2)*int(exp(0*1i*w/epsilo^2)*conj(p1{k}),w,0,s),[w,epsilo]);
end

for k=4:6
    m2{k}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*1i*w/epsilo^2)*(p1{k-3}- w*p{k-3}),w,0,s),[w,epsilo]);
    m2c{k}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*1i*w/epsilo^2)*conj(p1{k-3}- w*p{k-3}),w,0,s),[w,epsilo]);
    m0{k}=symfun(exp(-1i*s/epsilo^2)*int(exp(0*1i*w/epsilo^2)*(p1{k-3}- w*p{k-3}),w,0,s),[w,epsilo]);
    m0c{k}=symfun(exp(-1i*s/epsilo^2)*int(exp(0*1i*w/epsilo^2)*conj(p1{k-3}- w*p{k-3}),w,0,s),[w,epsilo]);
end
for k=1:6
    m2{k}=subs(m2{k},s,w);
    m2c{k}=subs(m2c{k},s,w);
    m0{k}=subs(m0{k},s,w);
    m0c{k}=subs(m0c{k},s,w);
end
for k=1:3
    r2{k}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*1i*w/epsilo^2)*(q2{k}-conj(q0c{k})),w,0,s),[w,epsilo]);
    r2c{k}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*1i*w/epsilo^2)*(conj(q2{k})-q0c{k}),w,0,s),[w,epsilo]);
    r0{k}=symfun(exp(-1i*s/epsilo^2)*int(exp(0*1i*w/epsilo^2)*(q2{k}-conj(q0c{k})),w,0,s),[w,epsilo]);
    r0c{k}=symfun(exp(-1i*s/epsilo^2)*int(exp(0*1i*w/epsilo^2)*(conj(q2{k})-q0c{k}),w,0,s),[w,epsilo]);
end
for k=4:6
    r2{k}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*1i*w/epsilo^2)*(q2c{k-3}-conj(q0{k-3})),w,0,s),[w,epsilo]);
    r2c{k}=symfun(exp(-1i*s/epsilo^2)*int(exp(2*1i*w/epsilo^2)*(conj(q2c{k-3})-q0{k-3}),w,0,s),[w,epsilo]);
    r0{k}=symfun(exp(-1i*s/epsilo^2)*int(exp(0*1i*w/epsilo^2)*(q2c{k-3}-conj(q0{k-3})),w,0,s),[w,epsilo]);
   r0c{k}=symfun(exp(-1i*s/epsilo^2)*int(exp(0*1i*w/epsilo^2)*(conj(q2c{k-3})-q0{k-3}),w,0,s),[w,epsilo]);
end
for k=1:6
    r2{k}=subs(r2{k},s,w);
    r2c{k}=subs(r2c{k},s,w);
    r0{k}=subs(r0{k},s,w);
    r0c{k}=subs(r0c{k},s,w);
end
savefile='symbol_intereal';
save(savefile,'p','p1','p2','q2','q2c','q0','q0c','q21','q2c1','q01','q0c1','r','m2','m2c','m0','m0c','r2','r2c','r0','r0c');














