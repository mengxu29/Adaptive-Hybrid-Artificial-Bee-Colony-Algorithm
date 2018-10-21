function y=qwblrnd(q,eta,beta,n)
% generate the random number with q weibull distribution
% use inverse transform sampling
if q == 1
    y = wblrnd(eta, beta, 1, n);
else if q < 2
        u = rand(1,n);
        y = eta*((1-(1-u).^((1-q)/(2-q)))/(1-q)).^(1/beta);
    else
        y = zeros(1,n);
    end
end

% % use rejection sampling method
% % http://en.wikipedia.org/wiki/Rejection_sampling#Adaptive_rejection_sampling
% % g(x): uniform dist
% % M: 1.1*max(f(x))
% 
% % n: the number of samples
% %% get sample range
% xmin = 1e-7;
% if (q<1)
%     xmax = yita*(1-q)^(-1/beta) * 0.1;
% else
%     xmax = yita*10;
% end
% %% get the max of pdf function
% M = max(pdf_qwbl(q,beta,yita,linspace(xmin, xmax, 1000))) * 1.1;
% %% generate
% nAccept = 1;
% samples = zeros(1, n);
% while (nAccept <= n)
%     x = xmin + (xmax - xmin) * rand;
%     u = rand;
%     if (u < pdf_qwbl(q,beta,yita,x)/M)
%         samples(nAccept) = x;
%         nAccept = nAccept + 1;
%     end
% end
% 
% y = samples;
% end
% 
% function y = pdf_qwbl(q,beta,yita,x)
% y = (2-q)*beta./yita*(x./yita).^(beta-1).*exp_q(-(x./yita).^beta,q);
% end
