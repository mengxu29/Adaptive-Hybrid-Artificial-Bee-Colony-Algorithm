B = 999;                                                                    
n = 21;
D = zeros(B+1,1);

LifeData=cell(1,1000);
load data_ApplicationExample_3_BCI_P_sample21;
LifeData{1}=[8,38,42,59,71,146,184,185,199,204,214,379,457,457,494,515,568,680,684,808,964];

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
p = (size(a,1) + 1) / (B + 1) ; 

