function y = Ie(x)
% % 
% % Euler = 0.5772156649;
% % EPS = 6E-08;
% % m = abs(log(EPS));
% % 
% % MAXIT = 100;
% % FPMIN = 1e-30;
% % 
% % ii=1:MAXIT;
% % if (x<FPMIN) 
% %     y=log(x)+Euler;
% % else if(x < -log(EPS))
% %         e=x.^ii./(ii.*factorial(ii));
% %         y = Euler + log(x) + sum(e);
% %     else
% %         e=factorial(ii)./x.^ii;
% %         y = (exp(x)/x)*(1+sum(e));
% %     end
% % end
% 
% if (~isreal(x))
%     errorFound = 1;
%     y=0;
%     return;
% end
if ~isreal(x)
    y=inf;
    return;
end

seg = exp([-150 -125 -100 -75 -55 -35 -21 -10]);
n_seg = length(seg);
I=find(x<seg);
if isempty(I)
    y=integral(@(t)exp(-t)./t,x,inf);
else if I(1) == 1%n_seg
        y=inf;%y=integral(@(t)exp(-t)./t,x,inf);
    else
        y=integral(@(t)exp(-t)./t,x,seg(I(1)));
        for ii = I(1):n_seg-1
            y=y+integral(@(t)exp(-t)./t,seg(ii),seg(ii+1));
        end
        y=y+integral(@(t)exp(-t)./t,seg(end),inf);
    end
end


%% use algorithm in textbook, by change the log into -log
% MAXIT = 100;
% EULER = 0.57721566;
% FPMIN = 1e-30;
% EPS = 6e-8;
% 
% if (x >= 0)
%     if (x < FPMIN) 
%         result = -log(x)+EULER;
%         if (result > 0.0)	
%             y= -log(x)+EULER;
%         else
%             y= FPMIN;
%         end
%     end
%     if (x <= -log(EPS)) %//test because of underflow.
%         sum=0.0; %//Use power series.
%         fact=1.0;
%         for k=1:MAXIT
%             fact = fact* x/k;
%             term=fact/k;
%             sum =sum+ term;
% %             if term<EPS*sum
% %                 break;
% %             end
%         end
% %         if (k >= MAXIT)
% %             printf('error in series');
% %         end%{ //nrerror("Series failed in ei");
%         result = sum-log(x)+EULER;
%         if (result > 0.0)
%             y= result;
%         else
%             y= FPMIN;
%         end
%     else %{ //Use asymptotic series.
%         sum=0.0; %//Start with second term.
%         term=1.0;
%         for k=1:MAXIT%(k=1;k<=MAXIT;k++) {
%             prev=term;
%             term = term* k/x;
%             if (term < EPS) 
%                 break;
%             end
%             %//Since final sum is greater than one, term itself approximates the relative error.
%             if (term < prev) 
%                 sum =sum+ term; %//Still converging: add new term.
%             else
%                 sum =sum- prev; %//Diverging: subtract previous term and
%                 break; %//exit.
%             end
%         end
%         exp_ = exp(x);
%         if (exp_ > 1e+303)
%             exp_ = 1e+303;
%         end
%         result = exp_*(1.0+sum)/x;
%         if (result > 0.0)
%             y= result;
%         else
%             y= FPMIN;
%         end
%     end
% else %%%%%%%%%%%%%%%%%%% never use this flow because the input is always positive
%     y = inf;
% end

%% use intergral for large
% MAXIT = 100;
% EULER = 0.57721566;
% FPMIN = 1e-30;
% EPS = 6e-8;
% DIVISION = 1e-3;
% 
% if (x < 0)
%     y=NAN;
% else if (x < FPMIN) 
%         result = - log(x)+EULER;
%         if (result > 0.0)	
%             y= -log(x)+EULER;
%         else
%             y= FPMIN;
%         end
%     else if (x <= DIVISION) %//test because of underflow.
%             sum=0.0; %//Use power series.
%             fact=1.0;
%             for k=1:MAXIT
%                 fact = fact* x/k;
%                 term=fact/k;
%                 sum =sum+ term;
%     %             if term<EPS*sum
%     %                 break;
%     %             end
%             end
%     %         if (k >= MAXIT)
%     %             printf('error in series');
%     %         end%{ //nrerror("Series failed in ei");
%             result = sum-log(x)+EULER;
%             if (result > 0.0)
%                 y= result;
%             else
%                 y= FPMIN;
%             end
%         else
%             y=integral(@(t)exp(-t)./t,x,inf);
%         end
%     end
% end

