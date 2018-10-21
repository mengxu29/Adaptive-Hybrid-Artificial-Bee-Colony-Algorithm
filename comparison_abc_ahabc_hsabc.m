load('LifeDataFile');

global lifes;

lifes = LifeData{18};

[Obj_AHABC, sol_AHABC, output_AHABC]=ABC_direct_throw_local_search_ahabc(1);
objMean_AHABC = mean(Obj_AHABC);
objStd_AHABC = std(Obj_AHABC);
solMean_AHABC = mean(sol_AHABC);
solStd_AHABC = std(sol_AHABC);

[Obj_ABC, sol_ABC, output_ABC]=ABC_direct_throw_local_search_ahabc(0);
objMean_ABC = mean(Obj_ABC);
objStd_ABC = std(Obj_ABC);
solMean_ABC = mean(sol_ABC);
solStd_ABC = std(sol_ABC);

[Obj_HSABC, sol_HSABC, output_HSABC]=ABC_direct_throw_local_search_hsabc(1);
objMean_HSABC = mean(Obj_HSABC);
objStd_HSABC = std(Obj_HSABC);
solMean_HSABC = mean(sol_HSABC);
solStd_HSABC = std(sol_HSABC);

save comparison_abc_ahabc_hsabc_18
