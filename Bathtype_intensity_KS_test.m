B = 999;                                                                    
n = 44;
D = zeros(B+1,1);

LifeData=cell(1,1000);
load data_Bathtype_intensity_BCI_P_qWeibull;

[TBF,~,~] = xlsread('bathtype_intensity.xlsx','sheet7','A1:A44');

lifes = zeros(1,length(TBF));
lifes(1) = TBF(1);
for i = 2:length(TBF)
    lifes(i) = lifes(i-1) + TBF(i) ;
end

LifeData{1} = lifes;

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

t = 0:0.01:20;
cdf = 1 - (1 - (1 - q).*(t./eta).^ beta).^((2 - q)/(1 - q)) ;
plot (t, cdf,'k--');
hold on;
t = LifeData{1};
cdfplot(t);
legend('Estimated CDF','Empirical CDF');
ylabel('Cumulative distribution function');
xlabel('Time (x1000 hours)');
title('');
axis([0 14 0 1]);
