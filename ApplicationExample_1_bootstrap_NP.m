global lifes;

LifeData = [0.058 0.070 0.090 0.105 0.113 0.121 0.153 0.159 0.224 0.421 0.570 0.596 0.618 0.834 1.019 1.104 1.497 2.027 2.234 2.372 2.433 2.505 2.690 2.877 2.879 3.166 3.455 3.551 4.378 4.872 5.085 5.272 5.341 8.952 9.188 11.399];
LifeData_resample=cell(1,1000);
n = 36;
solMean = [0.431837592891935 0.669675869909970 6.60868702451477];

for jj = 2:1000
    LifeData_resample{jj} = randsample(LifeData,n,true);
end

BCI_NP_sol = zeros(1,9);
BCI_NP_obj = zeros(1,3);

solutions = zeros(1000,3);
objection = zeros(1000,1);

solutions(1,:) = solMean;

for jj = 2:1000
    lifes = LifeData_resample{jj};
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
    
BCI_NP_sol(1) = sol_lb(1);
BCI_NP_sol(2) = sol_ub(1);
BCI_NP_sol(3) = sol_lb(2);
BCI_NP_sol(4) = sol_ub(2);
BCI_NP_sol(5) = sol_lb(3);
BCI_NP_sol(6) = sol_ub(3);
BCI_NP_sol(7) = sol_ub(1)-sol_lb(1);
BCI_NP_sol(8) = sol_ub(2)-sol_lb(2);
BCI_NP_sol(9) = sol_ub(3)-sol_lb(3);
BCI_NP_obj(1) = obj_lb;
BCI_NP_obj(2) = obj_ub;
BCI_NP_obj(3) = obj_ub - obj_lb;
    

save data_ApplicationExample_3_BCI_NP_test