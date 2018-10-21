load comparison_abc_ahabc_18;

% get the mean and std curve
t = 500:10:output{1}(end,1)-100;
nt = length(t);
objMeanCurve = zeros(1,nt);
objStdCurve = zeros(1,nt);
for kk = 1:length(t)
    obj_temp = zeros(1,length(output));
    for jj = 1:length(output)
        ind = find(output{jj}(:,1)<=t(kk),1,'last');   % find the cycle close to t(kk)
        obj_temp(jj) = output{jj}(ind,2);
    end
    objMeanCurve(kk) = mean(obj_temp);
    objStdCurve(kk) = std(obj_temp);
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

minObj =  min(min(Obj),min(Obj_ABC));
figure;

semilogy(t,(objMeanCurve-minObj),'k-',t,(objMeanCurve_ABC-minObj),'k--');

legend('AHABC','ABC');
title('Objective Convergence Speed');
xlabel('Number of function evaluations');
ylabel('The mean of max(lnL)-lnL in 30 runs');

% figure;
% plot(t,log(objStdCurve),'k-',t,log(objStdCurve_ABC),'r--');
% legend('AHABC','ABC');
% title('The Roubustness of Convergence');
% xlabel('Number of Function Evaluations');
% ylabel('Standard Deviation of Best Function Objective in 30 runs');


