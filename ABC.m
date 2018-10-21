function [GlobalMins,solutions]=ABC()
%/* Control Parameters of ABC algorithm*/
NP=200; %/* The number of colony size (employed bees+onlooker bees)*/
FoodNumber=NP/2; %/*The number of food sources equals the half of the colony size*/
limit=200; %/*A food source which could not be improved through "limit" trials is abandoned by its employed bee*/
maxCycle=5000; %/*The number of cycles for foraging {a stopping criteria}*/
maxBestTrial = 3000;
deltaObj = 1e-15;

global lifes;
%/* Problem specific variables*/

%objfun=@(x)qweibull_loglik(x); %cost function to be optimized
D=3; %/*The number of parameters of the problem to be optimized*/
max_eta = max(lifes)*2;
ub=[1.9 10 max_eta]; %/*lower bounds of the parameters. */
lb=[-10 0.00001 0.0001];%/*upper bound of the parameters.*/

runtime=10;%/*Algorithm can be run many times in order to see its robustness*/

%Foods [FoodNumber][D]; /*Foods is the population of food sources. Each row of Foods matrix is a vector holding D parameters to be optimized. The number of rows of Foods matrix equals to the FoodNumber*/
%ObjVal[FoodNumber];  /*f is a vector holding objective function values associated with food sources */
%Fitness[FoodNumber]; /*fitness is a vector holding fitness (quality) values associated with food sources*/
%trial[FoodNumber]; /*trial is a vector holding trial numbers through which solutions can not be improved*/
%prob[FoodNumber]; /*prob is a vector holding probabilities of food sources (solutions) to be chosen*/
%solution [D]; /*New solution (neighbour) produced by v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) j is a randomly chosen parameter and k is a randomly chosen solution different from i*/
%ObjValSol; /*Objective function value of new solution*/
%FitnessSol; /*Fitness value of new solution*/
%neighbour, param2change; /*param2change corresponds to j, neighbour corresponds to k in equation v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij})*/
%GlobalMin; /*Optimum solution obtained by ABC algorithm*/
%GlobalParams[D]; /*Parameters of the optimum solution*/
%GlobalMins[runtime]; /*GlobalMins holds the GlobalMin of each run in multiple runs*/

GlobalMins=zeros(1,runtime);
solutions = [];

for r=1:runtime
  
    fprintf('run time=%d \n',r);
%% initialize
% /*All food sources are initialized */
%/*Variables are initialized in the range [lb,ub]. If each parameter has different range, use arrays lb[j], ub[j] instead of lb and ub */

Range = repmat((ub-lb),[FoodNumber 1]);
Lower = repmat(lb, [FoodNumber 1]);
Foods = rand(FoodNumber,D) .* Range + Lower;
% initialze eta for all Foods
max_t = max(lifes);
% a = rand(FoodNumber,1)+1;
% mask = Foods(:,1)<1;
% Foods(mask,3) = (1-Foods(mask,1)).^(1./Foods(mask,2))*max_t.*a(mask);
% for ii = 1:FoodNumber
%     Foods(ii,:) = getNewFeasibleSolution(lb, ub, max_t);
% end

ObjVal = zeros(1,FoodNumber);
    for iFoods=1:FoodNumber
        ObjVal(iFoods)=objfun(Foods(iFoods,:));
    end
Fitness=calculateFitness(ObjVal);

%reset trial counters
trial=zeros(1,FoodNumber);

%/*The best food source is memorized*/
BestInd=find(ObjVal==min(ObjVal));
BestInd=BestInd(end);
GlobalMin=ObjVal(BestInd);
GlobalParams=Foods(BestInd,:);

iter=1;
bestTrial = 0;
preGlobalMin = GlobalMin;
while ((iter <= maxCycle)),
%% %%%%%%% EMPLOYED BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:(FoodNumber)
        
        %/*The parameter to be changed is determined randomly*/
        Param2Change=fix(rand*D)+1;
        
        %/*A randomly chosen solution is used in producing a mutant solution of the solution i*/
        neighbour=fix(rand*(FoodNumber))+1;
       
        %/*Randomly selected solution must be different from the solution i*/        
            while(neighbour==i)
                neighbour=fix(rand*(FoodNumber))+1;
            end;
        
       sol=Foods(i,:);
       %  /*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
       sol(Param2Change)=Foods(i,Param2Change)+(Foods(i,Param2Change)-Foods(neighbour,Param2Change))*(rand-0.5)*2;
        
       %  /*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
       sol = modifySolution(sol, lb, ub, max_t);
        
%         if (sol(1) < 1)
%             tmpEta = (1-sol(1))^(1/sol(2))*max_t;
%             if (tmpEta > sol(3))
% %                 sol(3) = min(tmpEta,ub(3));
%                   sol = getNewFeasibleSolution(lb, ub, max_t);
%             end
%         end
        
        %evaluate new solution
        ObjValSol=objfun(sol);
        FitnessSol=calculateFitness(ObjValSol);
        
       % /*a greedy selection is applied between the current solution i and its mutant*/
       PP = 0.1;%probability to choose bad solution
       if (isS1DominateS2(sol, Foods(i,:),FitnessSol,Fitness(i),lb,ub,max_t) || rand < PP) %/*If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i*/
            Foods(i,:)=sol;
            Fitness(i)=FitnessSol;
            ObjVal(i)=ObjValSol;
            trial(i)=0;
        else
            trial(i)=trial(i)+1; %/*if the solution i can not be improved, increase its trial counter*/
       end;
         
         
    end;
%%%%%%%%%%%%%%%%%%%%%%%% CalculateProbabilities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%/* A food source is chosen with the probability which is proportioal to its quality*/
%/*Different schemes can be used to calculate the probability values*/
%/*For example prob(i)=fitness(i)/sum(fitness)*/
%/*or in a way used in the method below prob(i)=a*fitness(i)/max(fitness)+b*/
%/*probability values are calculated by using fitness values and normalized by dividing maximum fitness value*/

prob=(0.9.*abs(Fitness)./max(abs(Fitness)))+0.1;
  
%% %%%%%%%%%%%%%%%%%%%%%% ONLOOKER BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=1;
t=0;
while(t<FoodNumber)
    if(rand<prob(i))
        t=t+1;
        %/*The parameter to be changed is determined randomly*/
        Param2Change=fix(rand*D)+1;
        
        %/*A randomly chosen solution is used in producing a mutant solution of the solution i*/
        neighbour=fix(rand*(FoodNumber))+1;
       
        %/*Randomly selected solution must be different from the solution i*/        
            while(neighbour==i)
                neighbour=fix(rand*(FoodNumber))+1;
            end;
        
       sol=Foods(i,:);
       %  /*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
       sol(Param2Change)=Foods(i,Param2Change)+(Foods(i,Param2Change)-Foods(neighbour,Param2Change))*(rand-0.5)*2;
        
       %  /*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
       sol = modifySolution(sol, lb, ub, max_t);
        
%         if (sol(1) < 1)
%             tmpEta = (1-sol(1))^(1/sol(2))*max_t;
%             if (tmpEta > sol(3))
% %                 sol(3) = min(tmpEta,ub(3));
%                   sol = getNewFeasibleSolution(lb, ub, max_t);
%             end
%         end
%         
        %evaluate new solution
        ObjValSol=objfun(sol);
        FitnessSol=calculateFitness(ObjValSol);
        
       % /*a greedy selection is applied between the current solution i and its mutant*/
       if (isS1DominateS2(sol, Foods(i,:),FitnessSol,Fitness(i),lb,ub,max_t)||rand<PP)%FitnessSol>Fitness(i) && isFeasible(sol, lb, ub, max_t)) %/*If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i*/
            Foods(i,:)=sol;
            Fitness(i)=FitnessSol;
            ObjVal(i)=ObjValSol;
            trial(i)=0;
        else
            trial(i)=trial(i)+1; %/*if the solution i can not be improved, increase its trial counter*/
       end;
    end;
    
    i=i+1;
    if (i==(FoodNumber)+1) 
        i=1;
    end;   
    
end; 

%/*The best food source is memorized*/
         indb=find(ObjVal==min(ObjVal));
         indb=indb(end);
         if (ObjVal(indb)<GlobalMin)
             GlobalMin=ObjVal(indb);
             GlobalParams=Foods(indb,:);
         end
         
         if (preGlobalMin > GlobalMin)
             bestTrial = 0; % reset bestTrial
         else 
             bestTrial = bestTrial + 1;
         end;
         
% stop critirion
         dObj = preGlobalMin - GlobalMin;
         if (((dObj ~= 0) && (dObj < deltaObj)) || (bestTrial > maxBestTrial))
             break;
         else
             preGlobalMin = GlobalMin;
         end
         
         
%% %%%%%%%%%% SCOUT BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%/*determine the food sources whose trial counter exceeds the "limit" value. 
%In Basic ABC, only one scout is allowed to occur in each cycle*/

ind=find(trial==max(trial));
ind=ind(end);
if (trial(ind)>limit)
    trial(ind)=0;
    sol=(ub-lb).*rand(1,D)+lb;
    a = rand()+1;
    if (sol(1) < 1)
%                 sol(3) = min(tmpEta,ub(3));
        sol = getNewFeasibleSolution(lb, ub, max_t); % initialize the eta value
    end
    ObjValSol=objfun(sol);
    FitnessSol=calculateFitness(ObjValSol);
    Foods(ind,:)=sol;
    Fitness(ind)=FitnessSol;
    ObjVal(ind)=ObjValSol;
    fprintf('scout happens at food %d\n',ind);
    %fprintf('best happens at food %d\n',indb);
end;

fprintf('iter=%d ObjVal=%g\n',iter,GlobalMin);
iter=iter+1;
end % End of ABC
GlobalMins(r)=GlobalMin;
solutions(end+1,:)=GlobalParams;
end; %end of runs
y = solutions;
save all

function fFitness=calculateFitness(fObjV)
fFitness=zeros(size(fObjV));
ind=find(fObjV>=0);
fFitness(ind)=1./(fObjV(ind)+1);
ind=find(fObjV<0);
fFitness(ind)=1+abs(fObjV(ind));

function sol = getNewFeasibleSolution(lb, ub, max_t)
D = 3;
a = rand(1,1)+1;
Foods = rand(1,D) .* (ub-lb) + lb;
while(Foods(1) < 1 && Foods(1,3) < (1-Foods(1,1)).^(1./Foods(1,2))*max_t)
    Foods = rand(1,D) .* (ub-lb) + lb;
end
sol = Foods;

function b = isFeasible(sol, lb, ub, max_t)
i1 = sum(sol < lb) + sum(sol > ub);
if i1 > 0
    b = 0;
else if ( sol(1) >= 1)
        b = 1;
    else if (sol(1,3) < (1-sol(1,1)).^(1./sol(1,2))*max_t)
            b = 0;
        else
            b = 1;
        end
    end
end

function b = isS1DominateS2(s1, s2,fit1, fit2, lb, ub, max_t)
i1 = sum(s1 < lb) + sum(s1 > ub);
i2 = sum(s2 < lb) + sum(s2 > ub);
if (s1(1) < 1)
    dq = (1-s1(1))^(1/s1(2))*max_t - s1(3);
    if (dq > 0)
        i1 = i1 + dq;
    end
end

if (s2(1) < 1)
    dq = (1-s2(1))^(1/s2(2))*max_t - s2(3);
    if (dq > 0)
        i2 = i2+dq;
    end
end

if (max(i1,i2) == 0)
    b = fit1 > fit2;
else
    b = i1<i2;
end

function new = modifySolution(sol, lb, ub, max_t)
ind=find(sol<lb);
sol(ind)=lb(ind);
ind=find(sol>ub);
sol(ind)=ub(ind);
% q_boundary = 1-(sol(3)./max_t).^sol(2);
% sol(1) = max(q_boundary, sol(1));
new = sol;

