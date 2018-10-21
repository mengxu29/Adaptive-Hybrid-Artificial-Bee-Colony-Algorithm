global lifes;

I = zeros(3,3);
ACI = zeros(1,9);

% lifes = [0.014 0.034 0.059 0.061 0.069 0.080 0.123 0.142 0.165 0.210 0.381 0.464 0.479 0.556 0.574 0.839 0.917 0.969 0.991 1.064 1.088 1.091 1.174 1.270 1.275 1.355 1.397 1.477 1.578 1.649 1.702 1.893 1.932 2.001 2.161 2.292 2.326 2.337 2.628 2.785 2.811 2.886 2.993 3.122 3.248 3.715 3.790 3.857 3.912 4.100 4.106 4.116 4.315 4.510 4.584 5.267 5.299 5.583 6.065 9.701];

lifes = [0.058 0.070 0.090 0.105 0.113 0.121 0.153 0.159 0.224 0.421 0.570 0.596 0.618 0.834 1.019 1.104 1.497 2.027 2.234 2.372 2.433 2.505 2.690 2.877 2.879 3.166 3.455 3.551 4.378 4.872 5.085 5.272 5.341 8.952 9.188 11.399];

[obj, sol] = ABC_direct_throw_local_search_adaptive(1);

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

% Z1 = -1.645;    % 90%
% Z2 = 1.645;

% Z1 = -1.282;    % 80%
% Z2 = 1.282;

Z1 = -1.036;    % 70% 
Z2 = 1.036;

    
q_ub = sol(1) + Z2*(var_q)^(0.5);
q_lb = sol(1) + Z1*(var_q)^(0.5);
beta_ub = sol(2) + Z2*(var_beta)^(0.5);
beta_lb = sol(2) + Z1*(var_beta)^(0.5);
eta_ub = sol(3) + Z2*(var_eta)^(0.5);
eta_lb = sol(3) + Z1*(var_eta)^(0.5);
    
ACI(1) = q_lb;
ACI(2) = q_ub;
ACI(3) = beta_lb;
ACI(4) = beta_ub;
ACI(5) = eta_lb;
ACI(6) = eta_ub;
ACI(7) = q_ub-q_lb;
ACI(8) = beta_ub-beta_lb;
ACI(9) = eta_ub-eta_lb;

save data_ApplicationExample_2_asymptotic_confidence_intervals_test
