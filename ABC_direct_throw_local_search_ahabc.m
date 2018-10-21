function [GlobalMins,solutions, plotOutPut]=ABC_direct_throw_local_search_ahabc(isLocalAdded)
%/* Control Parameters of ABC algorithm*/

NP=100; %/* The number of colony size (employed bees+onlooker bees)*/
FoodNumber=NP/2; %/*The number of food sources equals the half of the colony size*/
limit=150; %/*A food source which could not be improved through "limit" trials is abandoned by its employed bee*/
maxCycle=50000; %/*The number of cycles for foraging {a stopping criteria}*/
maxBestTrial = 1000;
maxnFcn = 5e5;
deltaObj = 1e-16;
Nscout = 0;

global nFcn;
global lifes;

D=3; %/*The number of parameters of the problem to be optimized*/
max_eta = mean(lifes);
ub=[1.9 10 max_eta]; %/*lower bounds of the parameters. */
lb=[-10 0.00001 0.0001];%/*upper bound of the parameters.*/

runtime = 5;%/*Algorithm can be run many times in order to see its robustness*/

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
       sol(Param2Change)=Foods(i,Param2Change)+(Foods(i,Param2Change)-Foods(neighbour,Param2Change))*(rand-0.5)*2;
        
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

%% local search

if (isLocalAdded)
    FitnessGlobalMin = calculateFitness(GlobalMin);
    fprintf('Nscout = %d ', Nscout);
    C = 1;
    Ns = C*Nscout*limit;

    [~,Indexs] = sort(ObjVal,'ascend');
    v = Foods(Indexs(1:D+1),:)';
    fv = ObjVal(Indexs(1:D+1));
    
    rho = 1; chi = 2; psi = 0.5; sigma = 0.5;
    onesn = ones(1,D);
    two2np1 = 2:D+1;
    one2n = 1:D;
    func_evals = 0;
    tolf = 1e-10;
    tolx = 1e-5;
    x = zeros(D,1);
    while func_evals < Ns
        if max(abs(fv(1)-fv(two2np1))) <= max(tolf,10*eps(fv(1))) && ...
                max(max(abs(v(:,two2np1)-v(:,onesn)))) <= max(tolx,10*eps(max(v(:,1))))
            break
        end

        % Compute the reflection point

        % xbar = average of the n (NOT n+1) best points
        xbar = sum(v(:,one2n), 2)/D;
        xr = (1 + rho)*xbar - rho*v(:,end);
        x(:) = xr; fxr = objfun(x);
        func_evals = func_evals+1;

        if fxr < fv(:,1)
            % Calculate the expansion point
            xe = (1 + rho*chi)*xbar - rho*chi*v(:,end);
            x(:) = xe; fxe = objfun(x);
            func_evals = func_evals+1;
            if fxe < fxr
                v(:,end) = xe;
                fv(:,end) = fxe;
                how = 'expand';
            else
                v(:,end) = xr;
                fv(:,end) = fxr;
                how = 'reflect';
            end
        else % fv(:,1) <= fxr
            if fxr < fv(:,D)
                v(:,end) = xr;
                fv(:,end) = fxr;
                how = 'reflect';
            else % fxr >= fv(:,n)
                % Perform contraction
                if fxr < fv(:,end)
                    % Perform an outside contraction
                    xc = (1 + psi*rho)*xbar - psi*rho*v(:,end);
                    x(:) = xc; fxc = objfun(x);
                    func_evals = func_evals+1;

                    if fxc <= fxr
                        v(:,end) = xc;
                        fv(:,end) = fxc;
                        how = 'contract outside';
                    else
                        % perform a shrink
                        how = 'shrink';
                    end
                else
                    % Perform an inside contraction
                    xcc = (1-psi)*xbar + psi*v(:,end);
                    x(:) = xcc; fxcc = objfun(x);
                    func_evals = func_evals+1;

                    if fxcc < fv(:,end)
                        v(:,end) = xcc;
                        fv(:,end) = fxcc;
                        how = 'contract inside';
                    else
                        % perform a shrink
                        how = 'shrink';
                    end
                end
                if strcmp(how,'shrink')
                    for j=two2np1
                        v(:,j)=v(:,1)+sigma*(v(:,j) - v(:,1));
                        x(:) = v(:,j); fv(:,j) = objfun(x);
                    end
                    func_evals = func_evals + D;
                end
            end
        end
        [fv,j] = sort(fv);
        v = v(:,j);
    end   % while
    
    for kk = 1:D+1
        if isFeasible(v(:,kk))
            Foods(Indexs(kk),:) = v(:,kk)';
            ObjVal(Indexs(kk)) = fv(kk);
        end
    end
end


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



