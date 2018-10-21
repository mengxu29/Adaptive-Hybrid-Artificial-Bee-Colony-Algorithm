q = 0.5;
beta = 0.5;
yita = 5;

t = [0:100];
y = (2-q)*(beta/yita)*(t./yita).^(beta-1).*exp_q(-(t./yita).^beta,q);

plot(t,y,'-');
xlabel('time');
ylabel('pdf');
title('f(t)');
legend('q=1.5 beta=0.5');