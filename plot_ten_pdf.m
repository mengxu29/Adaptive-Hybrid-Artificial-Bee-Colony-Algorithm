figure;
hold on

for jj = 1:10
    q = sol(jj,1);
    beta = sol(jj,2);
    yita = sol(jj,3);
        
    t = 1:yita/((1-q)^(1/beta));
    y = (2-q)*(beta/yita)*(t./yita).^(beta-1).*exp_q(-(t./yita).^beta,q);
%     R = (1-(1-q)*(t./yita).^beta).^((2-q)/(1-q));
%     h = y./R;
    plot(t,y,'b');
end

tmax = 420;
t=1:tmax;
beta=0.7125;
y=beta*t.^(beta-1)/tmax^beta;
plot(t,y,'r');
xlabel('time');
ylabel('pdf');
title('f(t)');

q = sol(:,1);
beta = sol(:,2);
yita = sol(:,3);
        
t_test = yita./((1-q).^(1./beta));

 




