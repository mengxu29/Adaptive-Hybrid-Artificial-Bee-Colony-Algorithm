figure;
xlabel('Time (x1000 hours)');
ylabel('');


% Reliability of Weibull
load data_Example2_Weibull;
eta = solMean(1);
beta = solMean(2);

t = 0:0.001:12; 
R_weibull = exp(-(t/eta).^beta);
plot(t,R_weibull,'r-');
hold;

% Conditional reliability of Weibull
m = 4;
t = m:0.001:12; 
Rm_weibull = exp(-(m/eta)^beta);
Rt_weibull = exp(-((t)/eta).^beta);
c_weibull = Rt_weibull/Rm_weibull;
plot(t,c_weibull,'r:');


m = 6;
t = m:0.001:12; 
Rm_weibull = exp(-(m/eta)^beta);
Rt_weibull = exp(-((t)/eta).^beta);
c_weibull = Rt_weibull/Rm_weibull;
plot(t,c_weibull,'r-.');

m = 8;
t = m:0.001:12; 
Rm_weibull = exp(-(m/eta)^beta);
Rt_weibull = exp(-((t)/eta).^beta);
c_weibull = Rt_weibull/Rm_weibull;
plot(t,c_weibull,'r--');
% q-Weibull
load comparison_ApplicationExample2;
q = solMean(1);
beta = solMean(2);
eta = solMean(3);
% Reliability of q-Weibull
t = 0:0.001:12; 
R_qweibull = (1-(1-q)*(t./eta).^beta).^((2-q)/(1-q));
plot(t,R_qweibull,'b-');
% Conditional reliability of q-Weibull
m = 4;
t = m:0.001:12; 
Rm_qweibull = (1-(1-q)*(m/eta).^beta).^((2-q)/(1-q));
Rt_qweibull = (1-(1-q)*((t)./eta).^beta).^((2-q)/(1-q));
c_qweibull = Rt_qweibull/Rm_qweibull;
plot(t,c_qweibull,'b:');

m = 6;
t = m:0.001:12; 
Rm_qweibull = (1-(1-q)*(m/eta).^beta).^((2-q)/(1-q));
Rt_qweibull = (1-(1-q)*((t)./eta).^beta).^((2-q)/(1-q));
c_qweibull = Rt_qweibull/Rm_qweibull;
plot(t,c_qweibull,'b-.');

m = 8;
t = m:0.001:12; 
Rm_qweibull = (1-(1-q)*(m/eta).^beta).^((2-q)/(1-q));
Rt_qweibull = (1-(1-q)*((t)./eta).^beta).^((2-q)/(1-q));
c_qweibull = Rt_qweibull/Rm_qweibull;
plot(t,c_qweibull,'b--');



xlabel('t (x1000 hours)');
ylabel('');
title('Conditional reliability');
legend('Weibull Reliability','Weibull Conditional Reliability with T=4','Weibull Conditional Reliability with T=6','Weibull Conditional Reliability with T=8','qWeibull Reliability','qWeibull Conditional Reliability with T=4','qWeibull Conditional Reliability with T=6','qWeibull Conditional Reliability with T=8');

axis([0 12.5 0 1.2]);
