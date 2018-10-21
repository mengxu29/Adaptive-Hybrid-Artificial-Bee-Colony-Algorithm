load('LifeDataFile');

global lifes;

lifes = LifeData{18};

[Obj, sol, output]=ABC_direct_throw_local_search_adaptive(1);
objMean = mean(Obj);
objStd = std(Obj);
solMean = mean(sol);
solStd = std(sol);

[Obj_ABC, sol_ABC, output_ABC]=ABC_direct_throw_local_search_adaptive(0);
objMean_ABC = mean(Obj_ABC);
objStd_ABC = std(Obj_ABC);
solMean_ABC = mean(sol_ABC);
solStd_ABC = std(sol_ABC);

save comparison_abc_ahabc_2
