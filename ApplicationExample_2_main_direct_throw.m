global lifes;

lifes =[0.478 0.583 0.753 0.753 0.801 0.834 0.944 0.959 1.377 1.534 2.400 2.639 2.944 2.981 3.392 3.393 3.904 4.829 5.328 5.562 6.122 6.331 6.531 11.019 12.986];

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
