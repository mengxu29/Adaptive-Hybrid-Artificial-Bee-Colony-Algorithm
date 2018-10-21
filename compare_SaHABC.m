load('LifeDataFile');

global lifes;

lifes = LifeData{18};

[Obj_SaHABC, sol_SaHABC, output_SaHABC]=ABC_direct_throw_SaHABC();
objMean_SaHABC = mean(Obj_SaHABC);
objStd_SaHABC = std(Obj_SaHABC);
solMean_SaHABC = mean(sol_SaHABC);
solStd_SaHABC = std(sol_SaHABC);
    
save compare_SaHABC_18
