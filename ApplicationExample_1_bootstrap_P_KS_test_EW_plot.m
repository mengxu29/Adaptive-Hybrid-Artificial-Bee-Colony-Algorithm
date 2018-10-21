alpha=sol(1);
theta=sol(2);
sigma=sol(3);
h=((alpha*theta/sigma)*((1-exp(-(t/sigma).^alpha)).^(theta-1)).*exp(-(t/sigma).^alpha).*(t/sigma).^(alpha-1))./(1-((1-exp(-(t/sigma).^alpha)).^theta))