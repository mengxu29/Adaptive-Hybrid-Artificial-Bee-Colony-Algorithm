global lifes;
% display('sampling life data...');

ES=[0.5 0.5 5 20;
    0.5 0.5 5 100;
    0.5 0.5 5 500;
    0.5 0.5 5 1000;
    1.5 0.5 5 20;
    1.5 0.5 5 100;
    1.5 0.5 5 500;
    1.5 0.5 5 1000;
    1 1 5 20;
    1 1 5 100;
    1 1 5 500;
    1 1 5 1000;
    0.5 1.5 5 20;
    0.5 1.5 5 100;
    0.5 1.5 5 500;
    0.5 1.5 5 1000;
    1.5 1.5 5 20;
    1.5 1.5 5 100;
    1.5 1.5 5 500;
    1.5 1.5 5 1000;];
objMean = zeros(20,1);
objStd = zeros(20,1);
solMean = zeros(20,3);
solStd = zeros(20,3);
for ii = 2    
    q = ES(ii,1);
    beta = ES(ii,2);
    eta = ES(ii,3);
    n = ES(ii,4);
    lifes = qwblrnd(q,eta,beta,n);
    [Obj, sol]=ABC();
    objMean(ii) = mean(Obj);
    objStd(ii) = std(Obj);
    solMean(ii,:) = mean(sol);
    solStd(ii,:) = std(sol);
end
save data_table3
%draw_profile_of_mle;