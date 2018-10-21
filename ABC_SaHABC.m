function [GlobalMins,solutions, plotOutPut]=ABC_SaHABC
%/* Control Parameters of ABC algorithm*/
NP=100; %/* The number of colony size (employed bees+onlooker bees)*/
FoodNumber=NP/2; %/*The number of food sources equals the half of the colony size*/
limit=150; %/*A food source which could not be improved through "limit" trials is abandoned by its employed bee*/
maxCycle=5000; %/*The number of cycles for foraging {a stopping criteria}*/
maxBestTrial = 1000;
maxnFcn = 2e5;
deltaObj = 1e-16;

global nFcn;
global lifes;

D=3; %/*The number of parameters of the problem to be optimized*/
max_eta = mean(lifes);
ub=[1.9 10 max_eta]; %/*lower bounds of the parameters. */
lb=[-10 0.00001 0.0001];%/*upper bound of the parameters.*/

runtime = 10;%/*Algorithm can be run many times in order to see its robustness*/

plotOutPut = cell(1,runtime);


GlobalMins=zeros(1,runtime);
solutions = [];

for r=1:runtime
  
    fprintf('run time=%d \n',r);
%% initialize
max_t = max(lifes);
nFcn = 0;
Nscout = 0;

Foods = zeros(FoodNumber,3);
  for i = 1:FoodNumber
      Foods(i,:) = getNewFeasibleSolution(lb, ub, max_t);
  end

ObjVal = zeros(1,FoodNumber);
    for iFoods=1:FoodNumber
        ObjVal(iFoods)=objfun(Foods(iFoods,:));
    end
Fitness=calculateFitness(ObjVal);

trial=zeros(1,FoodNumber);

%/*The best food source is memorized*/
BestInd=find(ObjVal==min(ObjVal));
BestInd=BestInd(end);
GlobalMin=ObjVal(BestInd);
GlobalParams=Foods(BestInd,:);

iter=1;
bestTrial = 0;
preGlobalMin = GlobalMin;
while (nFcn <= maxnFcn)%(iter <= maxCycle)),
    if (rand <= 0.5)
        fai = -exp(-3*iter/(25*maxCycle));
    else
        fai = exp(-3*iter/(25*maxCycle));
    end
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
       sol(Param2Change)=GlobalParams(Param2Change)+(Foods(i,Param2Change)-Foods(neighbour,Param2Change))*fai;
               
       if (~isFeasible(sol))
           continue;
       end

        ObjValSol=objfun(sol);
        FitnessSol=calculateFitness(ObjValSol);
               
        if (FitnessSol > Fitness(i))
            Foods(i,:)=sol;
            Fitness(i)=FitnessSol;
            ObjVal(i)=ObjValSol;
            trial(i)=0;
        else
            trial(i)=trial(i)+1; %/*if the solution i can not be improved, increase its trial counter*/
       end;
                 
    end;
%%%%%%%%%%%%%%%%%%%%%%%% CalculateProbabilities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prob=(0.9.*abs(Fitness)./max(abs(Fitness)))+0.1;
  
%% %%%%%%%%%%%%%%%%%%%%%% ONLOOKER BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=1;
t=0;
while(t<FoodNumber)
    if(rand<prob(i))
        t=t+1;
        
        Param2Change=fix(rand*D)+1;
        neighbour=fix(rand*(FoodNumber))+1;
        
            while(neighbour==i)
                neighbour=fix(rand*(FoodNumber))+1;
            end;
        
       sol=Foods(i,:);
       sol(Param2Change)=Foods(i,Param2Change)+(Foods(i,Param2Change)-Foods(neighbour,Param2Change))*(rand-0.5)*2+(GlobalParams(Param2Change)-Foods(i,Param2Change))*rand*1.5;
         
       if (~isFeasible(sol))
           continue;
       end
       
        ObjValSol=objfun(sol);
        FitnessSol=calculateFitness(ObjValSol);
       
       if (FitnessSol > Fitness(i))
            Foods(i,:)=sol;
            Fitness(i)=FitnessSol;
            ObjVal(i)=ObjValSol;
            trial(i)=0;
        else
            trial(i)=trial(i)+1;
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
%          dObj = preGlobalMin - GlobalMin;
%          if (((dObj ~= 0) && (dObj < deltaObj)) || (bestTrial > maxBestTrial))
%              fprintf('best trial is : %d \n',bestTrial);
%              break;c
%          else
%              preGlobalMin = GlobalMin;
%          end
                  
%% %%%%%%%%%% SCOUT BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind=find(trial==max(trial));
ind=ind(end);
if (trial(ind)>limit)
    trial(ind)=0;
    
    sol = getNewFeasibleSolution(lb, ub, max_t);
    Nscout = Nscout + 1;
    
    ObjValSol=objfun(sol);
    FitnessSol=calculateFitness(ObjValSol);
    Foods(ind,:)=sol;
    Fitness(ind)=FitnessSol;
    ObjVal(ind)=ObjValSol;
    fprintf('scout happens at food %d\n',ind);
end;


fprintf('iter=%d ObjVal=%g\n',iter,GlobalMin);
iter=iter+1;

plotOutPut{r}(end+1,:) = [nFcn, GlobalMin, mean(ObjVal), GlobalParams, Nscout]';
end % End of ABC


GlobalMins(r)=GlobalMin;
solutions(end+1,:)=GlobalParams;
end; %end of runs
save all

function fFitness=calculateFitness(fObjV)
fFitness=zeros(size(fObjV));
ind=find(fObjV>=0);
fFitness(ind)=1./(fObjV(ind)+1);
ind=find(fObjV<0);
fFitness(ind)=1+abs(fObjV(ind));

function b = isFeasible(sol)
global lifes;
if ( sol(1) < 2 && sol(2) > 0 && sol(3) > 0 && min((1-(1-sol(1))*(lifes/sol(3)).^sol(2)))>0 )
    b = 1;
else
    b = 0;
end


function sol = getNewFeasibleSolution(lb, ub, max_t)
D = 3;
Foods = rand(1,D).* (ub-lb) + lb;
while(Foods(1) < 1 && Foods(1,3) < (1-Foods(1,1)).^(1./Foods(1,2))*max_t)
    Foods = rand(1,D).* (ub-lb) + lb;
end
sol = Foods;



