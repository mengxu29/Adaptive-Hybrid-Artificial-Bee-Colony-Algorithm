B = 999;                                                                    
n = 60;
D = zeros(B+1,1);

LifeData=cell(1,1000);
load data_ApplicationExample_1_sample60_BCI_P;
LifeData{1}=[0.014 0.034 0.059 0.061 0.069 0.080 0.123 0.142 0.165 0.210 0.381 0.464 0.479 0.556 0.574 0.839 0.917 0.969 0.991 1.064 1.088 1.091 1.174 1.270 1.275 1.355 1.397 1.477 1.578 1.649 1.702 1.893 1.932 2.001 2.161 2.292 2.326 2.337 2.628 2.785 2.811 2.886 2.993 3.122 3.248 3.715 3.790 3.857 3.912 4.100 4.106 4.116 4.315 4.510 4.584 5.267 5.299 5.583 6.065 9.701];

Fn = zeros(1,n+1);
for j = 0:n
    Fn(j+1,1) = j / n ;
end
t = LifeData{1};
cdf = 1 - (1 - (1 - q).*(t./eta).^ beta).^((2 - q)/(1 - q)) ;
cdf = sort(cdf);
Dmais = max(abs(Fn(2:n+1) - cdf));
Dmenos = max(abs(Fn(1:n) - cdf));
D(1) = max(Dmais, Dmenos);

for jj = 2:1000
    
    t = LifeData{jj};
    cdf = 1 - (1 - (1 - solutions(jj,1)).*(t./solutions(jj,3)).^ solutions(jj,2)).^((2 - solutions(jj,1))/(1 - solutions(jj,1)));
    cdf = sort(cdf);
    Dmais = max(abs(Fn(2:n+1) - cdf));
    Dmenos = max(abs(Fn(1:n) - cdf));
    D(jj) = max(Dmais, Dmenos);
    
end

a = find(D > D(1)) ;
p = (size(a,1) + 1) / (B + 1) 

figure;

t = 0:0.001:12.4;
cdf = 1 - (1 - (1 - q).*(t./eta).^ beta).^((2 - q)/(1 - q)) ;
plot (t, cdf,'k--');
hold on;
t = LifeData{1};
cdfplot(t);
legend('Estimated CDF','Empirical CDF');
ylabel('Cumulative distribution function');
xlabel('Cycles (x1000)');
title('');
