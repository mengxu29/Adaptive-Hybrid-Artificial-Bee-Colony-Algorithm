% load data_AHSABC_AHABC;
% output_AHABC = output;
% Obj_AHABC = Obj;
% 
% load data_AHSABC_HSABC;
% output_HSABC = output;
% Obj_HSABC = Obj;
% 
% load data_AHSABC_ABC;
% output_ABC = output;
% Obj_ABC = Obj;

% get the mean and std curve
% t = 500:100:output_AHABC{1}(end,1)-100;
t = 500:100:200000;
nt = length(t);
objMeanCurve = zeros(1,nt);
objStdCurve = zeros(1,nt);
for kk = 1:length(t)
    obj_temp = zeros(1,length(output_AHABC));
    for jj = 1:length(output_AHABC)
        ind = find(output_AHABC{jj}(:,1)<=t(kk),1,'last');   % find the cycle close to t(kk)
        obj_temp(jj) = output_AHABC{jj}(ind,2);
    end
    objMeanCurve(kk) = mean(obj_temp);
    objStdCurve(kk) = std(obj_temp);
end

objMeanCurve_HSABC = zeros(1,nt);
objStdCurve_HSABC = zeros(1,nt);
for kk = 1:length(t)
    obj_temp = zeros(1,length(output_HSABC));
    for jj = 1:length(output_HSABC)
        ind = find(output_HSABC{jj}(:,1)<=t(kk),1,'last');   % find the cycle close to t(kk)
        obj_temp(jj) = output_HSABC{jj}(ind,2);
    end
    objMeanCurve_HSABC(kk) = mean(obj_temp);
    objStdCurve_HSABC(kk) = std(obj_temp);
end

objMeanCurve_ABC = zeros(1,nt);
objStdCurve_ABC = zeros(1,nt);
for kk = 1:length(t)
    obj_temp = zeros(1,length(output_ABC));
    for jj = 1:length(output_ABC)
        ind = find(output_ABC{jj}(:,1)<=t(kk),1,'last');   % find the cycle close to t(kk)
        obj_temp(jj) = output_ABC{jj}(ind,2);
    end
    objMeanCurve_ABC(kk) = mean(obj_temp);
    objStdCurve_ABC(kk) = std(obj_temp);
end


minObj =  min(min(min(Obj_AHABC),min(Obj_HSABC)),min(Obj_ABC));
figure;

semilogy(t,(objMeanCurve-minObj),'k-',t,(objMeanCurve_HSABC-minObj),'k--',t,(objMeanCurve_ABC-minObj),'r-');

legend('AHABC','HSABC','ABC');
title('Objective Convergence speed');
xlabel('Number of Function Evaluations');
ylabel('The mean of max(lnL)-lnL in 30 runs');

% figure;
% plot(t,log(objStdCurve),'k-',t,log(objStdCurve_ABC),'r--');
% legend('AHABC','ABC');
% title('The Roubustness of Convergence');
% xlabel('Number of Function Evaluations');
% ylabel('Standard Deviation of Best Function Objective in 30 runs');


