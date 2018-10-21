q =  1;
beta =  1.1458;
yita =  179.6559;

t0 = 0; % to simplify the problem, t0 is set to 0. But if you want you can add

global lifes;
t = lifes;
loglik_vector = log((2-q)*beta/(yita-t0)) + (beta-1)*log((t-t0)/(yita-t0)) + log(exp_q(-((t-t0)./(yita-t0)).^beta, q));
loglik = sum(loglik_vector)