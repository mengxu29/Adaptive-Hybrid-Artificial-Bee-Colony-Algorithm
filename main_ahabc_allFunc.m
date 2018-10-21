function [Obj, sol, output] = main_ahabc_allFunc()

Func_Paras.D = 4; %/*The number of parameters of the problem to be optimized*/
Func_Paras.ub = [200 10 1 1];%[1.9 10 max_eta]; %/*lower bounds of the parameters. */
Func_Paras.lb = [0 0 0 0];%[-10 0.00001 0.0001];%/*upper bound of the parameters.*/
Func_Paras.getNewFeasibleSolution = @()getNewFeasibleSolution(Func_Paras.lb,Func_Paras.ub,Func_Paras.D);
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
y = obj_IPRA(x);

function b = isFeasible(x)                              
if (x(1)>0 && x(2)>0 && x(3)>0 && x(4)>0 && x(3)<1 && x(4)<1)
    b = 1;
else
    b = 0;
end


function sol = getNewFeasibleSolution(lb,ub,D)
sol = rand(1,D).*(ub-lb) + lb;
