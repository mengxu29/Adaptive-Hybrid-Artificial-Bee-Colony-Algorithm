function loglik = qweibull_loglik(paras)
%% this function is to calculate the log of likelyhood, with a given paras
% paras: the parameters,i.e. [q, beta, yita]
% lifes: experimental test life data
q = paras(1);
beta = paras(2);
yita = paras(3);

t0 = 0; % to simplify the problem, t0 is set to 0. But if you want you can add

global lifes;
t = lifes;
loglik_vector = log((2-q)*beta/(yita-t0)) + (beta-1)*log((t-t0)/(yita-t0)) + log(exp_q(-((t-t0)./(yita-t0)).^beta, q));
loglik = sum(loglik_vector);
end
