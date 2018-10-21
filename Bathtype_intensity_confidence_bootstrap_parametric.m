global lifes;

LifeData=cell(1,1000);

q = 0.9495;
beta = 0.5652;
eta = 15.1495;
n = 44;
for jj = 2:1000
    LifeData{jj} = qwblrnd(q,eta,beta,n);
end

BCI_P_sol = zeros(1,9);
BCI_P_obj = zeros(1,3);

solutions = zeros(1000,3);
objection = zeros(1000,1);
    
solutions(1,:) = [q beta eta];

for jj = 2:1000
    lifes = LifeData{jj};
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
    
BCI_P_sol(1) = sol_lb(1);
BCI_P_sol(2) = sol_ub(1);
BCI_P_sol(3) = sol_lb(2);
BCI_P_sol(4) = sol_ub(2);
BCI_P_sol(5) = sol_lb(3);
BCI_P_sol(6) = sol_ub(3);
BCI_P_sol(7) = sol_ub(1)-sol_lb(1);
BCI_P_sol(8) = sol_ub(2)-sol_lb(2);
BCI_P_sol(9) = sol_ub(3)-sol_lb(3);
BCI_P_obj(1) = obj_lb;
BCI_P_obj(2) = obj_ub;
BCI_P_obj(3) = obj_ub-obj_lb;

save data_Bathtype_intensity_BCI_P_qWeibull
