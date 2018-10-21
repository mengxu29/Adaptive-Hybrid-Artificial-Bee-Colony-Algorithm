global lifes;
% %runMe;
% display('using fminsearch to optimize');
% fminsearch(@(x)objfun((x)),[1;0.5;2.5])
load('LifeDataFile');

objMean = zeros(20,1);
objStd = zeros(20,1);
solMean = zeros(20,3);
solStd = zeros(20,3);
for ii = 1  
    q = ES(ii,1);
    beta = ES(ii,2);
    eta = ES(ii,3);
    n = ES(ii,4);
    lifes = LifeData{ii};
    tic;
    [Obj, sol]=ABC_modify_eta();
    toc;
    objMean(ii) = mean(Obj);
    objStd(ii) = std(Obj);
    solMean(ii,:) = mean(sol);
    solStd(ii,:) = std(sol);
end
save data_ABC_modify_eta_1
%draw_profile_of_mle;