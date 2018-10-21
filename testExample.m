q = -60481.8096601520;
beta = 0.712403241159020;
eta = 2318514349.83386;
n=200;

LifeData = qwblrnd(q,eta,beta,n);

global lifes;
lifes = LifeData;

[Obj, sol]=ABC_direct_throw_local_search();
objMean = mean(Obj);
objStd = std(Obj);
solMean = mean(sol);
solStd = std(sol);
    
save data_ApplicationExample_1_ABC_direct_throw