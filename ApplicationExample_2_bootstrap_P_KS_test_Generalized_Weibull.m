function ApplicationExample_2_bootstrap_P_KS_test_Generalized_Weibull()
global lifes;
lifes = [0.058 0.070 0.090 0.105 0.113 0.121 0.153 0.159 0.224 0.421 0.570 0.596 0.618 0.834 1.019 1.104 1.497 2.027 2.234 2.372 2.433 2.505 2.690 2.877 2.879 3.166 3.455 3.551 4.378 4.872 5.085 5.272 5.341 8.952 9.188 11.399];

%% solve the mle of original data

Func_Paras.D = 3; %/*The number of parameters of the problem to be optimized*/
Func_Paras.ub = [5 5 5];%[1.9 10 max_eta]; %/*lower bounds of the parameters. */
Func_Paras.lb = [0 0 0];%[-10 0.00001 0.0001];%/*upper bound of the parameters.*/
Func_Paras.getNewFeasibleSolution = @()getNewFeasibleSolution(Func_Paras.lb,Func_Paras.ub,Func_Paras.D);
Func_Paras.anyObjFun = @negative_likelihood_func;
Func_Paras.isFeasible = @isFeasible;

ABC_Paras.NP=100; %/* The number of colony size (employed bees+onlooker bees)*/
ABC_Paras.FoodNumber=ABC_Paras.NP/2; %/*The number of food sources equals the half of the colony size*/
ABC_Paras.limit=150; %/*A food source which could not be improved through "limit" trials is abandoned by its employed bee*/
ABC_Paras.maxCycle=50000; %/*The number of cycles for foraging {a stopping criteria}*/
ABC_Paras.maxBestTrial = 1000;
ABC_Paras.maxnFcn = 5e5;
ABC_Paras.deltaObj = 1e-16;
ABC_Paras.isLocalAdded = 1;
ABC_Paras.method = 'AHABC';
ABC_Paras.C = 1;
ABC_Paras.Ns = 5;
ABC_Paras.runtime = 5;%/*Algorithm can be run many times in order to see its robustness*/

[obj, sol, ~]=ABC_ahabc_allFunc(Func_Paras,ABC_Paras);
solMean = mean(sol);
solstd = std(sol);
objMean = mean(obj);
objstd = std(obj);

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

save data_ApplicationExample_1_bootstrap_P_KS_test_Generalized_weibull

%% define the MLE problems
function neg_loglik=negative_likelihood_func(paras)
global nFcn;
nFcn = nFcn+1;

global lifes;
t = lifes;

alpha = paras(1);
theta = paras(2);
lamda = paras(3);

neg_loglik_vector = -(log(alpha)+log(theta)+(theta-1)*log(t)+((1/lamda)-1)*log(1-alpha*lamda*(t.^theta)));
neg_loglik = sum(neg_loglik_vector);

function samples = randSampleFunc(paras,n)
alpha = paras(1);
theta = paras(2);
lamda = paras(3);
u = rand(1,n);
samples = ((1-(1-u).^lamda)/(alpha*lamda)).^(1/theta);

function b = isFeasible(paras)
if(paras(3)>0 && 11.399 > (paras(1)*paras(3))^(-1/paras(2)))
    b = 0;
else
    b = paras(1)>0 && paras(2)>0;
end


function sol = getNewFeasibleSolution(lb,ub,D)
sol = rand(1,D).*(ub-lb) + lb;
while ~isFeasible(sol)
    sol = rand(1,D).*(ub-lb) + lb;
end

