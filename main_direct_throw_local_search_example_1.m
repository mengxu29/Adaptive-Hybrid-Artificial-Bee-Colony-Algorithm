global lifes;

load('LifeDataFile');

lifes = [5 11 21 31 46 75 98 122 145 165 196 224 245 293 321 330 350 420];
[Obj, sol]=ABC_direct_throw();
objMean = mean(Obj);
objStd = std(Obj);
solMean = mean(sol);
solStd = std(sol);
    
save data_ApplicationExample_1_ABC_direct_throw