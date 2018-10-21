% explore the property of MLE function
global lifes;
lifes = [5 11 21 31 46 75 98 122 145 165 196 224 245 293 321 330 350 420];

[qq, ee]=meshgrid(-100:0.01:-99,2e5:1e2:2.75e5);
[n1, n2]=size(qq);
ff = zeros(n1,n2);
for jj=1:n1
    for kk = 1:n2
        ff(jj,kk) = -qweibull_loglik([qq(jj,kk) 0.7124 ee(jj,kk)]);
    end
end
figure;
mesh(qq,ee,ff);
xlabel('q');
ylabel('eta');
zlabel('f');