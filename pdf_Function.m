
q = [0.5 0.5 1 1.5 1.5];
beta = [0.5 1.5 1 0.5 1.5];
eta = [5 5 5 5 5];

t = 0:0.2:20;
y = zeros(5,length(t));
R = zeros(5,length(t));
h = zeros(5,length(t));
h(3,:) = 0.2*ones(1,length(t));


for ii = 1:5
    y(ii,:) = (2-q(ii))*(beta(ii)/eta(ii))*(t./eta(ii)).^(beta(ii)-1).*exp_q(-(t./eta(ii)).^beta(ii),q(ii));
end

for ii = 1:5
    R(ii,:)= (1-(1-q(ii))*(t./eta(ii)).^beta(ii)).^((2-q(ii))/(1-q(ii)));
end

for ii = 1:5
    h(ii,:) = y(ii,:)./R(ii,:);
end


h(h>2)=inf;
h2 = h(2,:); h2(h2==0)=inf; h(2,:) = h2;
h(3,:) = 0.2*ones(1,length(t));

plot(t,h,'-');
xlabel('t');
ylabel('h_q(t)');
title('');
legend('q=0.5 \beta=0.5','q=0.5 \beta=1.5','q=1 \beta=1','q=1.5 \beta=0.5','q=1.5 \beta=1.5');

axis([0 20 0 1]);