global lifes;
load('data_ABC_direct_throw_local_search');

LifeData=cell(20,1000);
for ii=1:20
    q = ES(ii,1);
    beta = ES(ii,2);
    eta = ES(ii,3);
    n = ES(ii,4);
    for jj = 2:1000
        LifeData{ii,jj} = qwblrnd(q,eta,beta,n);
    end
end

MSE = zeros(20,3);
bias = zeros(20,3);
variance = zeros(20,3);
solutions = zeros(1000,3);

for ii=8
    solutions(1,:) = solMean(ii,:);
    for jj = 2:1000
        lifes = LifeData{ii,jj};
        [Obj, sol] = ABC_direct_throw_local_search();
        solutions(jj,:) = sol;
        
        fprintf('sample number=%d \n',jj);
    end
    bias(ii,:) = mean(solutions)-ES(ii,1:3);
    variance(ii,:) = var(solutions);
    MSE(ii,:) = bias(ii,:).^2 + var(solutions);
end

save data_MSE_8_20000
