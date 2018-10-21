

% q-Weibull
load comparison_ApplicationExample2;
q = solMean(1);
beta = solMean(2);
eta = solMean(3);

B = 0.0318474160696793093233;
mean = eta*(2+(1/(1-q))+(1/beta))*((1-q)^(-1/beta))*B

