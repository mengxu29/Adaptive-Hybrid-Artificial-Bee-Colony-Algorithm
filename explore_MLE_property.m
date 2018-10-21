% explore the property of MLE function
global lifes;

load('LifeDataFile');

ii = 20;
lifes = LifeData{ii};

[qq, bb]=meshgrid(-10:0.005:1.9,0.1:0.05:10);
[n1, n2]=size(qq);
ff = zeros(n1,n2);
for jj=1:n1
    for kk = 1:n2
        eta = 5;
        ff(jj,kk) = -qweibull_loglik([qq(jj,kk) bb(jj,kk) eta]);
    end
end

mesh(qq,bb,ff);
xlabel('q');
ylabel('beta');
zlabel('f');