% q = 1.15214826302140;
% beta = 1.29510510206534;
% yita = 2.99519424753865;


 
q = 0.4318;
beta = 0.6697;
yita = 6.6087;


t = 0:0.001:15;
y = (2-q)*(beta/yita)*(t./yita).^(beta-1).*exp_q(-(t./yita).^beta,q);
R = (1-(1-q)*(t./yita).^beta).^((2-q)/(1-q));
h = y./R;

% figure;
% plot(t,y, 'k');
% xlabel('Cycles (x1000)');
% ylabel('Probability density function');
% 
% axis([0 8 0 1.6]);

% figure;
% plot(t,R, 'k');
% xlabel('Cycles (x1000)');
% ylabel('Reliability');
% axis([0 14 0 1.1]);
 
 
figure;
plot(t,h, 'k');
xlabel('Time (x1000 hours)');
ylabel('Hazard rate function');

axis([0 17 0 5]);