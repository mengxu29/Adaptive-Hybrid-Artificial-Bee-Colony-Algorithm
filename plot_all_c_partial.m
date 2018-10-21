load comparison_all_data;
figure;
semilogy(t,(objMeanCurve_ABC-minObj),...
t(1:length(objMeanCurve_C1)),(objMeanCurve_C1-minObj),...
t(1:length(objMeanCurve_C125)),(objMeanCurve_C125-minObj));

legend('C=0','C=1','C=125');
title('Effect of C on Convergence Speed');
xlabel('Number of function evaluations');
ylabel('The mean of max(lnL)-lnL in 30 runs');
xlim([0 200000]);

figure;
semilogy(t,(objStdCurve_ABC),...
t(1:length(objStdCurve_C1)),(objStdCurve_C1),...
t(1:length(objStdCurve_C125)),(objStdCurve_C125));

legend('C=0','C=1','C=125');
title('Effect of C on Convergence Variability');
xlabel('Number of function evaluations');
ylabel('The std of max(lnL)-lnL in 30 runs');
xlim([0 200000]);

