global lifes;
load('data_ABC_direct_throw_local_search');

LifeData_resample=cell(20,1000);
for ii=1:20
    n = ES(ii,4);
    for jj = 2:1000
        LifeData_resample{ii,jj} = randsample(LifeData{ii},n,true);
    end
end

BCI_NP_sol = zeros(20,9);
BCI_NP_obj = zeros(20,3);

solutions = zeros(1000,3);
objection = zeros(1000,1);
sol_lb = zeros(1,3);
sol_ub = zeros(1,3);

for ii=8
    
    solutions(1,:) = solMean(ii,:);
    objection(1) = objMean(ii);
    for jj = 2:1000
        lifes = LifeData_resample{ii,jj};
        [Obj, sol] = ABC_direct_throw_local_search();
        solutions(jj,:) = sol;  
        objection(jj) = Obj;
        
        fprintf('sample number=%d \n',jj);
    end
    y = sort(solutions);
    z = sort(objection);
    sol_lb = y(51,:);
    sol_ub = y(950,:);
    obj_lb = z(51);
    obj_ub = z(950);
    
    BCI_NP_sol(ii,1) = sol_lb(1);
    BCI_NP_sol(ii,2) = sol_ub(1);
    BCI_NP_sol(ii,3) = sol_lb(2);
    BCI_NP_sol(ii,4) = sol_ub(2);
    BCI_NP_sol(ii,5) = sol_lb(3);
    BCI_NP_sol(ii,6) = sol_ub(3);
    BCI_NP_sol(ii,7) = sol_ub(1) - sol_lb(1);
    BCI_NP_sol(ii,8) = sol_ub(2) - sol_lb(2);
    BCI_NP_sol(ii,9) = sol_ub(3) - sol_lb(3);
    BCI_NP_obj(ii,1) = obj_lb;
    BCI_NP_obj(ii,2) = obj_ub;
    BCI_NP_obj(ii,3) = obj_ub - obj_lb;
    
end

save data_bootstrap_confidence_intervals_NP_8