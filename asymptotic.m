global lifes;
load('LifeDataFile');

I = zeros(3,3);
var = zeros(3,3);
ACI = zeros(20,9);

for ii = 1:20
    lifes = LifeData{ii};
    [obj, sol] = ABC_direct_throw_local_search();
    
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
    
    q_beta_loglik_vector = log(t./eta)./((eta./t).^beta-1+q).^2;
    q_beta_loglik = sum(q_beta_loglik_vector);
    
    beta_eta_loglik_vector = (1/eta)*(-1+1./(((eta./t).^beta)-1+q))-beta*(eta^(beta-1))*(log(eta./t)./((t.^beta).*((eta./t).^beta-1+q).^2));
    beta_eta_loglik = sum(beta_eta_loglik_vector);
    
    eta_q_loglik_vector = -(beta/eta)*(((eta./t).^beta)-1+q).^(-2);
    eta_q_loglik = sum(eta_q_loglik_vector);
    
    I(1,1) = q_loglik;
    I(2,2) = beta_loglik;
    I(3,3) = eta_loglik;
    I(1,2) = q_beta_loglik;
    I(2,1) = q_beta_loglik;
    I(1,3) = eta_q_loglik;
    I(3,1) = eta_q_loglik;
    I(2,3) = beta_eta_loglik;
    I(3,2) = beta_eta_loglik;
    
    var = -I^(-1);
    var_q = var(1,1);
    var_beta = var(2,2);
    var_eta = var(3,3);
    
    Z1 = -1.64;
    Z2 = 1.64;
    
    q_ub = sol(1) + Z2*(var_q)^(0.5);
    q_lb = sol(1) + Z1*(var_q)^(0.5);
    beta_ub = sol(2) + Z2*(var_beta)^(0.5);
    beta_lb = sol(2) + Z1*(var_beta)^(0.5);
    eta_ub = sol(3) + Z2*(var_eta)^(0.5);
    eta_lb = sol(3) + Z1*(var_eta)^(0.5);
    
    ACI(ii,1) = q_lb;
    ACI(ii,2) = q_ub;
    ACI(ii,3) = beta_lb;
    ACI(ii,4) = beta_ub;
    ACI(ii,5) = eta_lb;
    ACI(ii,6) = eta_ub;
    ACI(ii,7) = q_ub-q_lb;
    ACI(ii,8) = beta_ub-beta_lb;
    ACI(ii,9) = eta_ub-eta_lb;
    
    
end
 
save data_asymptotic_confidence_intervals
