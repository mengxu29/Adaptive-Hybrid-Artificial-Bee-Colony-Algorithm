q = 0.5;
beta = 0.772598259573171;
yita = 4455.20187855435;


% t = yita/(1-q)^(1/beta);

t = 0:0.01:10920;
y = (2-q)*(beta/yita)*(t./yita).^(beta-1).*exp_q(-(t./yita).^beta,q);
R = (1-(1-q)*(t./yita).^beta).^((2-q)/(1-q));
h = y./R;

figure;
plot(t,y, 'k');
xlabel('Distance (x1000km)');
ylabel('Probability density function');

% axis([0 30 0 0.25]);

figure;
plot(t,R, 'k');
xlabel('Distance (x1000km)');
ylabel('Reliability');
% axis([0 30 0 1.1]);


figure;
plot(t,h, 'k');
xlabel('Distance (x1000km)');
ylabel('Hazard rate function');

% axis([0 30 0 0.05]);