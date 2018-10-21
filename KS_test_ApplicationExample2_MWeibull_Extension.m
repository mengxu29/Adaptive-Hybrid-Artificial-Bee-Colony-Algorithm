B = 999;                                                                    
D = zeros(B+1,1);

LifeData=cell(1,1000);
load data_ApplicationExample_2_BCI_P_MWeibull_Extension;

LifeData{1}=[0.058 0.070 0.090 0.105 0.113 0.121 0.153 0.159 0.224 0.421 0.570 0.596 0.618 0.834 1.019 1.104 1.497 2.027 2.234 2.372 2.433 2.505 2.690 2.877 2.879 3.166 3.455 3.551 4.378 4.872 5.085 5.272 5.341 8.952 9.188 11.399];

n = length(LifeData{1});
Fn = zeros(1,n+1);
for j = 0:n
    Fn(j+1,1) = j / n ;
end
t = LifeData{1};
cdf = 1-exp(solutions(1,1)*solutions(1,3)*(1-exp((t/solutions(1,1)).^solutions(1,2))));
cdf = sort(cdf);
Dmais = max(abs(Fn(2:n+1) - cdf));
Dmenos = max(abs(Fn(1:n) - cdf));
D(1) = max(Dmais, Dmenos);

for jj = 2:1000
    
    t = LifeData{jj};
    cdf = 1-exp(solutions(jj,1)*solutions(jj,3)*(1-exp((t/solutions(jj,1)).^solutions(jj,2))));
    cdf = sort(cdf);
    Dmais = max(abs(Fn(2:n+1) - cdf));
    Dmenos = max(abs(Fn(1:n) - cdf));
    D(jj) = max(Dmais, Dmenos);
    
end

a = find(D > D(1)) ;
p = (size(a,1) + 1) / (B + 1) 

figure;

t = 0:0.001:12.4;
cdf = 1-exp(solutions(1,1)*solutions(1,3)*(1-exp((t/solutions(1,1)).^solutions(1,2))));
plot (t, cdf,'k--');
hold on;
t = LifeData{1};
cdfplot(t);
legend('Estimated CDF','Empirical CDF');
ylabel('Cumulative distribution function');
xlabel('Time (x1000 hours)');
title('MWeibull Extension');
