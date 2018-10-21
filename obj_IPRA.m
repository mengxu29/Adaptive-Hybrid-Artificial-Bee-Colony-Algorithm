function obj = obj_IPRA(paras)
%% this function is to calculate the log of likelyhood, with a given paras

alpha = paras(1);
beta = paras(2);
q_x = paras(3);
q_z = paras(4);
r = 0.647;

% exclude qx qz < 0
if (q_x<0 || q_z < 0 || q_x > 1 || q_z > 1 || beta < 0 || alpha < 0 || r < 0 || r > 1)
    obj = inf;
    return
end

% example 1
% j = [1 0 0 1 0 0 1 0 1 1 1 1 0 1 0 1 1 1 1 1 1 0 1 1];
% y = [41 65 10 31 1 146 32 14 1 1 26 23 43 28 12 29 28 1 21 11 43 2 70 1];

% example 2
j = [1 1 1 1 1 1 1 0 1 0 1 1 1 0 0 0 1 1 1 1 1 1 0 0 0 1 1 1 1 1 1 1 1 0 1 0 0 0 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 1 1 0 1 1 1 0 0 1 1 0 1 1 1 1];
y = [220 13 1 6 25 5 3 6 6 2 7 1 5 25 3 5 32 3 1 12 36 1 11 10 4 1 1 32 14 1 12 7 28 10 24 8 1 1 1 19 2 1 1 13 6 3 6 2 12 1 3 7 2 12 12 117 3 4 2 2 30 97 65 47 7 18 8 80 61 11 28 12 13 24 3 10 4 85 28 5 76 49 4 32 17];

n=length(y);

u = (j*q_z+(1-j)*q_x).*y;
w = zeros(1,n);
loglik_vector = zeros(1,n);

w(1) = u(1);
for i = 2:n
    w(i) = w(i-1)+u(i);
end


loglik_vector(1) = n*log(beta)-n*beta*log(alpha)+(beta-1)*log(y(1))+log((1-r)*exp(-(y(1)/alpha)^beta)+r*Ie((y(1)/alpha)^beta));
for i = 2:n
    loglik_vector(i)=(beta-1)*log(y(i)+w(i-1))+log((1-r)*exp(-((y(i)+w(i-1))/alpha)^beta)+r*Ie(((y(i)+w(i-1))/alpha)^beta))-log(exp(-(w(i-1)/alpha)^beta)-r*((w(i-1)/alpha)^beta)*Ie((w(i-1)/alpha)^beta));
end
loglik = sum(loglik_vector);
obj = -loglik;
if(obj==-Inf)
    isFound=1;
end
