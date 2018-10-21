figure;
xlabel('Time(x1000 hours)');
ylabel('CDF');
t = 0:0.001:15.37;

% Weibull
load data_Example2_Weibull;
cdf = 1 - exp(-(t/solMean(1)).^solMean(2));
plot(t,cdf,'k-.','LineWidth',1.2);
hold;

% q-Weibull
load comparison_ApplicationExample2;
q = solMean(1);
beta = solMean(2);
eta = solMean(3);

R = (1-(1-q)*(t./eta).^beta).^((2-q)/(1-q));
cdf = 1-R;
plot(t,cdf,'k-','LineWidth',1.2);

% ENH Lemonte
load data_Example2_ENH;
cdf = (1 - exp(1-(1+(solMean(3)*t)).^solMean(1))).^solMean(2);
plot(t,cdf,'k--','LineWidth',1.2);

% Modified Weibull Extension
load data_Example2_MWeibull_Extension;
cdf = 1-exp(solMean(1)*solMean(3)*(1-exp((t/solMean(1)).^solMean(2))));
plot(t,cdf,'r:','LineWidth',1.4);

% Empirical CDF
LifeData{1} = [0.058 0.070 0.090 0.105 0.113 0.121 0.153 0.159 0.224 0.421 0.570 0.596 0.618 0.834 1.019 1.104 1.497 2.027 2.234 2.372 2.433 2.505 2.690 2.877 2.879 3.166 3.455 3.551 4.378 4.872 5.085 5.272 5.341 8.952 9.188 11.399];
t = LifeData{1};
n = 36;
Fn = zeros(1,n+1);
for j = 0:n
    Fn(j+1,1) = j / n ;
end
cdfplot(t);

xlabel('Time(x1000 hours)');
ylabel('CDF');
title('');
legend('Weibull','q-Weibull','ENH','Modified Weibull Extension','Empirical CDF');

axis([0 17 0 1]);


