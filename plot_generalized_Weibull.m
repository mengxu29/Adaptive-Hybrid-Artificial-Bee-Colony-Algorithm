h=alpha*theta*(t.^(theta-1)).*((1-alpha*lamda*(t.^theta)).^(-1))


 f = (2-q)*(beta/eta)*(t/eta).^(beta-1).*exp_q(-(t./eta).^beta,q);
 
 R= (1-(1-q)*(t/eta).^beta).^((2-q)/(1-q));
 h=f./R
 
 
 plot(t,h,'-');
xlabel('t');
ylabel('h_q(t)');
title('');
legend('Generalized Weibull','q-Weibull');

axis([0 20 0 1]);