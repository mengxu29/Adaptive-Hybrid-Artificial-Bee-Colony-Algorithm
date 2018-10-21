load comparison_abc_ahabc;

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
plot(t,log(objMeanCurve-minObj),'k-',t,log(objMeanCurve_ABC-minObj),'r--');
legend('AHABC','ABC');
title('Objective Convergency speed');
xlabel('Number of Function Evaluations');
ylabel('The mean value of Best Function Objective in 10 runs');

figure;
plot(t,log(objStdCurve),'k-',t,log(objStdCurve_ABC),'r--');
legend('AHABC','ABC');
title('The Roubustness of Convergency');
xlabel('Number of Function Evaluations');
ylabel('Standard Deviation of Best Function Objective in 10 runs');

