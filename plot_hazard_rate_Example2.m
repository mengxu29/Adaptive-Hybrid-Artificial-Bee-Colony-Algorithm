figure;
xlabel('Time(x1000 hours)');
ylabel('Hazard rate function');
t = 0:0.001:15.37;

% Weibull
load data_Example2_Weibull;
eta = solMean(1);
beta = solMean(2);

h = (beta/(eta^beta))*t.^(beta-1);

plot(t,h,'k-.','LineWidth',1.2);
hold;

% q-Weibull
load comparison_ApplicationExample2;
q = solMean(1);
beta = solMean(2);
eta = solMean(3);

y = (2-q)*(beta/eta)*(t./eta).^(beta-1).*exp_q(-(t./eta).^beta,q);
R = (1-(1-q)*(t./eta).^beta).^((2-q)/(1-q));
h = y./R;

plot(t,h,'k-','LineWidth',1.2);

% ENH Lemonte
load data_Example2_ENH;
alpha = solMean(1);
beta = solMean(2);
lamda = solMean(3);
a(1)=alpha;
a(2)=beta;
a(3)=lamda;

h=a(1)*a(2)*a(3).*(1+a(3)*t).^(a(1)-1).*exp(1-(1+a(3)*t).^a(1)).*(1-exp(1-(1+a(3)*t).^a(1))).^(a(2)-1)./(1-(1-exp(1-(1+a(3)*t).^a(1))).^a(2));
plot(t,h,'k--','LineWidth',1.2);

% Modified Weibull Extension
load data_Example2_MWeibull_Extension;
alpha = solMean(1);
beta = solMean(2);
lamda = solMean(3);

h = lamda*beta*(t/alpha).^(beta-1).*exp((t/alpha).^beta);
plot(t,h,'r:','LineWidth',1.4);



xlabel('Time(x1000 hours)');
ylabel('Hazard rate function');
title('');
legend('Weibull','q-Weibull','ENH','Modified Weibull Extension');

axis([0 17 0 4]);



%Empirical
% 
% t =[0.014 0.034 0.059 0.061 0.069 0.080 0.123 0.142 0.165 0.210 0.381 0.464 0.479 0.556 0.574 0.839 0.917 0.969 0.991 1.064 1.088 1.091 1.174 1.270 1.275 1.355 1.397 1.477 1.578 1.649 1.702 1.893 1.932 2.001 2.161 2.292 2.326 2.337 2.628 2.785 2.811 2.886 2.993 3.122 3.248 3.715 3.790 3.857 3.912 4.100 4.106 4.116 4.315 4.510 4.584 5.267 5.299 5.583 6.065 9.701];
% emp_h = zeros(1,60);
% emp_h(1)=1/(t(1)*61);
% 
% for i=2:60
%     emp_h(i)=1/((t(i)-t(i-1))*(62-i));
% end
% 
% plot(t,emp_h,'g-.');


% % Generalized Weibull
% load data_Example2_Generalized_Weibull;
% alpha = solMean(1);
% theta = solMean(2);
% lamda = solMean(3);
% 
% h = alpha*theta*(t.^(theta-1)).*((1-alpha*lamda*(t.^theta)).^(-1));
% plot(t,h,'r-');
