figure;
[X Y]=meshgrid(-10:0.1:1.9,0.1:0.02:1);
[m n]=size(X);
Z=zeros(m,n);
for i = 1:m*n
    Z(i) = qweibull_loglik([X(i),Y(i),2.5]);
end
surf(X,Y,Z)
xlabel('q')
ylabel('\beta')
title('Likelyhood with \eta = 2.5')

figure;
[X Y]=meshgrid(0.07:0.02:1.9,0.1:0.1:5);
[m n]=size(X);
Z=zeros(m,n);
for i = 1:m*n
    Z(i) = qweibull_loglik([0.9,X(i),Y(i)]);
end
surf(X,Y,Z)
xlabel('\beta')
ylabel('\eta')
title('Likelyhood with q = 0.9')

figure;
[X Y]=meshgrid(-10:0.1:1.5,0.1:0.1:5);
[m n]=size(X);
Z=zeros(m,n);
for i = 1:m*n
    Z(i) = qweibull_loglik([X(i),0.5,Y(i)]);
end
surf(X,Y,Z)
xlabel('q')
ylabel('\eta')
title('Likelyhood with \beta = 0.5')