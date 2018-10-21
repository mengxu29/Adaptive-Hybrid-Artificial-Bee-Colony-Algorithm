global lifes;
load('data_ABC_local_search');

LifeData=cell(20,1000);
q = solMean(2,1);
beta = solMean(2,2);
eta = solMean(2,3);
n = 20;
for jj = 2:1000
    LifeData{ii,jj} = qwblrnd(q,eta,beta,n);
end

BCI_P_sol = zeros(20,9);
BCI_P_obj = zeros(20,3);

solutions = zeros(1000,3);
objection = zeros(1000,1);
    
solutions(1,:) = solMean(2,:);
objection(1) = objMean(2);
for jj = 2:1000
    lifes = LifeData{ii,jj};
    [Obj, sol] = ABC_direct_throw_local_search_adaptive(1);
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

BCI_P_sol(ii,1) = sol_lb(1);
BCI_P_sol(ii,2) = sol_ub(1);
BCI_P_sol(ii,3) = sol_lb(2);
BCI_P_sol(ii,4) = sol_ub(2);
BCI_P_sol(ii,5) = sol_lb(3);
BCI_P_sol(ii,6) = sol_ub(3);
BCI_P_sol(ii,7) = sol_ub(1) - sol_lb(1);
BCI_P_sol(ii,8) = sol_ub(2) - sol_lb(2);
BCI_P_sol(ii,9) = sol_ub(3) - sol_lb(3);
BCI_P_obj(ii,1) = obj_lb;
BCI_P_obj(ii,2) = obj_ub;
BCI_P_obj(ii,3) = obj_ub - obj_lb;
    

save data_bootstrap_confidence_intervals_test

