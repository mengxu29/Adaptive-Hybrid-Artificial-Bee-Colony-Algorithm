function y = objfun(x)
global nFcn;
nFcn = nFcn + 1;
y = -qweibull_loglik(x);