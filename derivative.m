t = lifes;
q = sol(1);
beta = sol(2);
eta = sol(3);

q_loglik_vector = -(2-q)^(-2)+2*(1-q)^(-3).*log(1-(1-q).*((t./eta).^beta))+2*(1-q)^(-2)*(((eta./t).^beta)-1+q).^(-1)-(1/(1-q)).*((eta./t).^beta-1+q).^(-2);
q_loglik = sum(q_loglik_vector);

beta_loglik_vector = -1/(beta^2)+(eta^beta)*((log(t./eta).*log(eta./t))./((t.^beta).*(((eta./t).^beta)-1+q).^2));
beta_loglik = sum(beta_loglik_vector);

eta_loglik_vector = (beta/(eta^2)).*(1-(1./(((eta./t).^beta)-1+q)))-(beta^2)*(eta^(beta-2))./((t.^beta).*(((eta./t).^beta)-1+q).^2);
eta_loglik = sum(eta_loglik_vector);
  
    
    