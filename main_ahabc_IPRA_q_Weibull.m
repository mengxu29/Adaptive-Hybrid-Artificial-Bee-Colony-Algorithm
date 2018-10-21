function [Obj, sol, output] = main_ahabc_IPRA_q_Weibull()
global lifes;

% lifes = [41 65 10 31 1 146 32 14 1 1 26 23 43 28 12 29 28 1 21 11 43 2 70 1];  
lifes = [220 13 1 6 25 5 3 6 6 2 7 1 5 25 3 5 32 3 1 12 36 1 11 10 4 1 1 32 14 1 12 7 28 10 24 8 1 1 1 19 2 1 1 13 6 3 6 2 12 1 3 7 2 12 12 117 3 4 2 2 30 97 65 47 7 18 8 80 61 11 28 12 13 24 3 10 4 85 28 5 76 49 4 32 17];

max_eta = mean(lifes);
Func_Paras.max_t = max(lifes);
Func_Paras.D = 6; %/*The number of parameters of the problem to be optimized*/
Func_Paras.ub = [1.9 10 max_eta 1 1 1];%[1.9 10 max_eta]; %/*lower bounds of the parameters. */
Func_Paras.lb = [-10 0.00001 0.00001 0 0 0];%[-10 0.00001 0.0001];%/*upper bound of the parameters.*/
Func_Paras.getNewFeasibleSolution = @()getNewFeasibleSolution(Func_Paras.lb,Func_Paras.ub,Func_Paras.D,Func_Paras.max_t);
Func_Paras.anyObjFun = @anyObjFun;
Func_Paras.isFeasible = @isFeasible;

ABC_Paras.NP=100; %/* The number of colony size (employed bees+onlooker bees)*/
ABC_Paras.FoodNumber=ABC_Paras.NP/2; %/*The number of food sources equals the half of the colony size*/
ABC_Paras.limit=150; %/*A food source which could not be improved through "limit" trials is abandoned by its employed bee*/
ABC_Paras.maxCycle=1000; %/*The number of cycles for foraging {a stopping criteria}*/
ABC_Paras.maxBestTrial = 1000;
ABC_Paras.maxnFcn = 1e5;
ABC_Paras.deltaObj = 1e-16;
ABC_Paras.isLocalAdded = 1;
ABC_Paras.method = 'AHABC';
ABC_Paras.C = 1;
ABC_Paras.Ns = 5;
ABC_Paras.runtime = 1;%/*Algorithm can be run many times in order to see its robustness*/

tic;
[Obj, sol, output]=ABC_ahabc_allFunc(Func_Paras,ABC_Paras);
toc;
objMean = mean(Obj);
objStd = std(Obj);
solMean = mean(sol);
solStd = std(sol);
dataFileName = ['data_IPRA' ABC_Paras.method];
save(dataFileName)


function y=anyObjFun(x)
global nFcn;
nFcn = nFcn+1;
y = obj_IPRA_q_Weibull(x);


function b = isFeasible(x)
global lifes;
if ( x(1) < 2 && x(2) > 0 && x(3) > 0 && min((1-(1-x(1))*(lifes/x(3)).^x(2)))>0 && x(4)>0 && x(5)>0 && x(6)>0 && x(4)<1 && x(5)<1 && x(6)<1)
    b = 1;
else
    b = 0;
end


function sol = getNewFeasibleSolution(lb, ub, D, max_t)
x = rand(1,D).* (ub-lb) + lb;
while(x(1) < 1 && x(1,3) < (1-x(1,1)).^(1./x(1,2))*max_t)
    x = rand(1,D).* (ub-lb) + lb;
end
sol = x;

