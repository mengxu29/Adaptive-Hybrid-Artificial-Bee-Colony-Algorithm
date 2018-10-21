
global lifes;

% lifes = [0.014 0.034 0.059 0.061 0.069 0.080 0.123 0.142 0.165 0.210 0.381 0.464 0.479 0.556 0.574 0.839 0.917 0.969 0.991 1.064 1.088 1.091 1.174 1.270 1.275 1.355 1.397 1.477 1.578 1.649 1.702 1.893 1.932 2.001 2.161 2.292 2.326 2.337 2.628 2.785 2.811 2.886 2.993 3.122 3.248 3.715 3.790 3.857 3.912 4.100 4.106 4.116 4.315 4.510 4.584 5.267 5.299 5.583 6.065 9.701];
lifes = [0.058 0.070 0.090 0.105 0.113 0.121 0.153 0.159 0.224 0.421 0.570 0.596 0.618 0.834 1.019 1.104 1.497 2.027 2.234 2.372 2.433 2.505 2.690 2.877 2.879 3.166 3.455 3.551 4.378 4.872 5.085 5.272 5.341 8.952 9.188 11.399];
% lifes = [0.01 0.01 0.01 0.01 0.01 0.01 0.02 0.02 0.02 0.02 0.03 0.04 0.06 0.08 0.10 0.10 0.12 0.12 0.12 0.13 0.14 0.15 0.15 0.15 0.16 0.16 0.17 0.18 0.18 0.19 0.20 0.21 0.22 0.23 0.25 0.26 0.28 0.28 0.30 0.32 0.34 0.36 0.38 0.39 0.41 0.41 0.42 0.43 0.44 0.44 0.45 0.45 0.50 0.53 0.56 0.58 0.58 0.61 0.62 0.62 0.62 0.64 0.66 0.70 0.70 0.70 0.72 0.77 0.78 0.78 0.80 0.82 0.83 0.85 0.86 0.96 0.97 0.98 0.99 1.05 1.06 1.07 1.18 1.35 1.36 1.42 1.55 1.59 1.65 1.73 1.77 1.79 1.80 1.91 2.09 2.14 2.15 2.15 2.31 2.33 2.36 2.36 2.43 2.45 2.50 2.51 2.58 2.64 2.68 3.08 3.94 4.12 4.33 4.42 4.53 4.88 4.97 5.11 5.32 5.55 6.63 6.89 7.62 11.41 11.76 11.85 12.36 13.22];

[Obj, sol, output]=ABC_direct_throw_local_search_adaptive(1);
objMean = mean(Obj);
objStd = std(Obj);
solMean = mean(sol);
solStd = std(sol);


save comparison_Example;

q = solMean(1);
beta = solMean(2);
eta = solMean(3);
t = 0:0.001:12.4;
y = (2-q)*(beta/eta)*(t./eta).^(beta-1).*exp_q(-(t./eta).^beta,q);
R = (1-(1-q)*(t./eta).^beta).^((2-q)/(1-q));
h = y./R;

plot(t,h,'-');

xlabel('t');
ylabel('h_q(t)');
title('');

axis([0 14 0 5]);