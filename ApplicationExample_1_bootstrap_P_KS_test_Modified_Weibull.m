function ApplicationExample_1_bootstrap_P_KS_test_Modified_Weibull()
global lifes;
lifes = [0.014 0.034 0.059 0.061 0.069 0.080 0.123 0.142 0.165 0.210 0.381 0.464 0.479 0.556 0.574 0.839 0.917 0.969 0.991 1.064 1.088 1.091 1.174 1.270 1.275 1.355 1.397 1.477 1.578 1.649 1.702 1.893 1.932 2.001 2.161 2.292 2.326 2.337 2.628 2.785 2.811 2.886 2.993 3.122 3.248 3.715 3.790 3.857 3.912 4.100 4.106 4.116 4.315 4.510 4.584 5.267 5.299 5.583 6.065 9.701];

%% solve the mle of original data

Func_Paras.D = 3; %/*The number of parameters of the problem to be optimized*/
Func_Paras.ub = [10 10 10];%[1.9 10 max_eta]; %/*lower bounds of the parameters. */
Func_Paras.lb = [0 0 0];%[-10 0.00001 0.0001];%/*upper bound of the parameters.*/
Func_Paras.getNewFeasibleSolution = @()getNewFeasibleSolution(Func_Paras.lb,Func_Paras.ub,Func_Paras.D);
Func_Paras.anyObjFun = @negative_likelihood_func;
Func_Paras.isFeasible = @isFeasible;

ABC_Paras.NP=100; %/* The number of colony size (employed bees+onlooker bees)*/
ABC_Paras.FoodNumber=ABC_Paras.NP/2; %/*The number of food sources equals the half of the colony size*/
ABC_Paras.limit=150; %/*A food source which could not be improved through "limit" trials is abandoned by its employed bee*/
ABC_Paras.maxCycle=50000; %/*The number of cycles for foraging {a stopping criteria}*/
ABC_Paras.maxBestTrial = 1000;
ABC_Paras.maxnFcn = 2e5;
ABC_Paras.deltaObj = 1e-16;
ABC_Paras.isLocalAdded = 1;
ABC_Paras.method = 'AHABC';
ABC_Paras.C = 1;
ABC_Paras.Ns = 5;
ABC_Paras.runtime = 10;%/*Algorithm can be run many times in order to see its robustness*/

[obje, sol, ~]=ABC_ahabc_allFunc(Func_Paras,ABC_Paras);
solMean = mean(sol);

%% sample the rest 999 samples
LifeData=cell(1,1000);
LifeData{1} = lifes;

n = length(lifes);
for jj = 2:1000
    % define the sample function for different pdf
    LifeData{jj} = randSampleFunc(solMean,n);
end

solutions = zeros(1000,Func_Paras.D);
objection = zeros(1000,1);
    
solutions(1,:) = solMean;

%% solve MLE for these 999 samples
for jj = 2:1000
    lifes = LifeData{jj};
    
    ABC_Paras.runtime = 1;% only run once for one sampe, and return one solution
    [Obj, sol, ~]=ABC_ahabc_allFunc(Func_Paras,ABC_Paras);
    solutions(jj,:) = sol;
    objection(jj) = Obj;
         
    fprintf('sample number=%d \n',jj);
end

save data_ApplicationExample_1_bootstrap_P_KS_test_Modified_Weibull

%% define the MLE problems
function neg_loglik=negative_likelihood_func(paras)
global nFcn;
nFcn = nFcn+1;

global lifes;
t = lifes;

a = paras(1);
b = paras(2);
lamda = paras(3);

neg_loglik_vector = -(log(a)+log(b+lamda*t)+(b-1)*log(t)+lamda*t-a*t.^b.*exp(lamda*t));
neg_loglik = sum(neg_loglik_vector);

function samples = randSampleFunc(paras,n)
eta = paras(1);
beta = paras(2);
samples = qwblrnd(1,eta,beta,n);

function b = isFeasible(paras)
b = paras(1)>0 && paras(2)>0 && paras(3)>0;

function sol = getNewFeasibleSolution(lb,ub,D)
sol = rand(1,D).*(ub-lb) + lb;
while ~isFeasible(sol)
    sol = rand(1,D).*(ub-lb) + lb;
end

