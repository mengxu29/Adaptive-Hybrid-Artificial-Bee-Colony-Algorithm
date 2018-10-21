q = 1.157;
beta = 1.299;
yita = 2.969;

t0 = 0; % to simplify the problem, t0 is set to 0. But if you want you can add

t = [0.478 0.583 0.753 0.753 0.801 0.834 0.944 0.959 1.377 1.534 2.4 2.639 2.944 2.981 3.392 3.393 3.904 4.829 5.328 5.562 6.122 6.331 6.531 11.019 12.986];
loglik_vector = log((2-q)*beta/(yita-t0)) + (beta-1)*log((t-t0)/(yita-t0)) + log(exp_q(-((t-t0)./(yita-t0)).^beta, q));
loglik = sum(loglik_vector)