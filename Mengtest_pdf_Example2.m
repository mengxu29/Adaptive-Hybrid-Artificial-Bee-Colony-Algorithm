q = 1.15214826302140;
beta = 1.29510510206534;
yita = 2.99519424753865;

t = 0:50;
y = (2-q)*(beta/yita)*(t./yita).^(beta-1).*exp_q(-(t./yita).^beta,q);
R = (1-(1-q)*(t./yita).^beta).^((2-q)/(1-q));
h = y./R;

figure;
plot(t,y);
xlabel('time(days)');
ylabel('probability density function');


figure;
plot(t,R);
xlabel('time(days)');
ylabel('Reliability');


figure;
plot(t,h);
xlabel('time(days)');
ylabel('failure rate');
