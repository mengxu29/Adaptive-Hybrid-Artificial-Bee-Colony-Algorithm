global lifes;
lifes = [5 11 21 31 46 75 98 122 145 165 196 224 245 293 321 330 350 420];
t=  lifes;
beta= 0.712;
tmax=max(t);
n=length(t);
N = 10000;
xx=linspace(-0.021,0.02,N);
yy=zeros(1,N);
for jj = 1:N
    yy(jj)=n^2-(n+sum(log(1-xx(jj)*t.^beta)))*sum(1./(1-xx(jj)*t.^beta));
end

figure;
plot(xx,yy);

x_true=fzero(@(x) (n^2-(n+sum(log(1-x*t.^beta)))*sum(1./(1-x*t.^beta))), -0.0205)

s=sum(log(1-x_true*t.^beta));
q=1+s/(n+s)
eta=((1-q)/x_true)^(1/beta)
xlimit=(1/x_true)^(1/beta)
figure;
ezplot(@(x,beta) (n^2-(n+sum(log(1-x*t.^beta)))*sum(1./(1-x*t.^beta))));

figure;
kk = linspace(0.01,2,N);
yy=zeros(1,N);
for jj=1:N
    k = kk(jj);
    yy(jj)=sum(t.^k.*log(t))/sum(t.^k)-1/n*sum(log(t))-1/k;
end
plot(kk,yy);

beta=fzero(@(k)sum(t.^k.*log(t))/sum(t.^k)-1/n*sum(log(t))-1/k,1.15)
eta=(sum(t.^beta)/n)^(1/beta)