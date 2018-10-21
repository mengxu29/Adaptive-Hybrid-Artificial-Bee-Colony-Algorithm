function [y]=hjorthrand(n,delta,beta,theta)
% n: the number of samples
% delta, beta, theta: parameters 

% refer to online page:
% https://www.jstor.org/stable/1268388?seq=1#page_scan_tab_contents

%% define the cdf function for random sampling
% notice that this is a general method to sampling any distribution
% numerically. 
% if you want to try other cdf, just change the following lines.
cdf_function = @(t) 1 - exp(-delta*t^2/2)/((1+beta*t)^(theta/beta));

% check the parameter for hjorth distribution
if (theta <=0 || delta <= 0 || beta <= 0)
    y = INF;
    return;
end

%% use fminsearch to solve the inverse function
y = zeros(n,1);
for ii = 1:n
    u = rand();
    options = optimset('TolFun',1e-15,'MaxFunEvals',2000);
    x0 = 0.1; % the initial guess of the t
    [y(ii),fval] = fminsearch(@(t) abs(cdf_function(t)-u), x0, options);
    if (fval>1e-10)
        display('something goes wrong with the sampling function, please check');
    end
end

