global lifes;
% display('sampling life data...');
% beta = 0.5;
% eta = 5;
% q = 1.5;
% n = 1000;
% %lifes=qwblrnd(q, beta, yita, n);
% lifes = qwblrnd(q,eta,beta,n);
% hist(lifes)
% %runMe;
% display('using fminsearch to optimize');
% fminsearch(@(x)objfun((x)),[1;0.5;2.5])
LifeData=cell(1,20);
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
for ii = 1:20  
    q = ES(ii,1);
    beta = ES(ii,2);
    eta = ES(ii,3);
    n = ES(ii,4);
    LifeData{ii} = qwblrnd(q,eta,beta,n);
end
%save LifeDataFile
%draw_profile_of_mle;